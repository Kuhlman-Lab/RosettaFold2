"""
Methods for generating single and paired MSAs for RosettaFold2.
"""

import os
import requests
import time
import random
import tarfile
import copy
import numpy as np
import logging
import shutil
import sys
sys.path.append('/proj/kuhl_lab/alphafold/')
from alphafold.data import parsers
from alphafold.data import pipeline
from alphafold.data import pipeline_multimer
from alphafold.data import feature_processing
from alphafold.data import msa_pairing
from alphafold.data.tools import hhsearch
from alphafold.common import protein
from alphafold.notebooks import notebook_utils

from typing import Sequence, Optional, Dict, Tuple, MutableMapping, Union

sys.path.append('/proj/kuhl_lab/alphafold/run/')
from utils import utils

logger = logging.getLogger('features')

# (filename, [sequence])
CleanQuery = Tuple[str, str]

# {sequence: (raw MSA, raw templates)}
RawInput = Dict[str, Tuple[str, str]]


def getRawInputs(
        queries: Sequence[CleanQuery],
        msa_mode: str,
        use_filter: bool = True,
        use_templates: bool = False,
        custom_msa_path: Optional[str] = None,
        insert_msa_gaps: bool = False,
        update_msa_query_seq: float = 1.00,
        output_dir: str = '',
        use_pairing: bool = False,
        design_run: bool = False,
        proc_id: Optional[int] = None) -> RawInput:
    """ Computes and gathers raw a3m lines for the list of queries. """
    
    raw_inputs = {}
    
    # Gather unique sequences to run MMseqs2 in a batch.
    unique_sequences = []
    for query in queries:
        filename = query[0]
        seqs = query[1]
        
        for seq in seqs:
            if seq not in unique_sequences:
                unique_sequences.append(seq)
                
    # If a custom msas provided, parse and keep track of them.
    custom_msas = {}
    if custom_msa_path is not None:
        custom_msas.update(getCustomMSADict(custom_msa_path))


        new_custom_msas = {}
        for msa_seq in custom_msas:
            for input_seq in unique_sequences:
 
                if msa_seq in input_seq and insert_msa_gaps:
                    n_gaps = len(input_seq) - len(msa_seq)
                    if n_gaps > 0:
                        prepend_append = input_seq.split(msa_seq)
                        assert len(prepend_append) == 2, "Uh.. something weird at the insert gap spot."
                        n_prepend = len(prepend_append[0])
                        n_append = len(prepend_append[1])
                        old_a3m = custom_msas[msa_seq]
                        old_a3m_lines = old_a3m.split('\n')

                        new_lines = update_msa_lines(old_a3m_lines, input_seq, n_prepend, n_append)
                        new_custom_msas[input_seq] = '\n'.join(new_lines)

                elif len(msa_seq) == len(input_seq) and update_msa_query_seq:
                    seq_sim = sum([s1 == s2 for s1, s2 in zip(list(msa_seq), list(input_seq))]) / len(msa_seq)

                    if seq_sim > update_msa_query_seq:
                        old_a3m = custom_msas[msa_seq]
                        old_a3m_lines = old_a3m.split('\n')

                        new_lines = update_msa_lines(old_a3m_lines, input_seq, 0, 0)
                        new_custom_msas[input_seq] = '\n'.join(new_lines)
                        
        custom_msas.update(new_custom_msas)

        # Make sure our custom msas only includes sequences that were inputted
        final_custom_msas = {}
        for seq in unique_sequences:
            if seq in custom_msas:
                final_custom_msas[seq] = custom_msas[seq]
        custom_msas = final_custom_msas
            

    # If not using templates and custom MSA provided, remove sequence from
    # MMseqs2 queue.
    if not use_templates and custom_msas != {}:
        for sequence in custom_msas:
            if sequence in unique_sequences:
                unique_sequences.remove(sequence)
    print(unique_sequences)
    # If we want MSAs or need templates, let's run MMseqs2.            
    if msa_mode != 'single_sequence' and unique_sequences != []:
        use_env = True if msa_mode == 'MMseqs2-U+E' else False

        if proc_id is not None:
            prefix = f'{os.path.join(output_dir, "mmseqs2")}_{proc_id}'
        else:
            prefix = f'{os.path.join(output_dir, "mmseqs2")}'
        #print('features::getRawInputs:', prefix)
        a3m_lines, template_paths = runMMseqs2(
            prefix=prefix,
            sequences=unique_sequences,
            use_env=use_env,
            use_filter=use_filter,
            use_templates=use_templates,
            use_pairing=use_pairing)
        if design_run:
            if use_filter:
                mode = 'env' if use_env else 'all'
            else:
                mode = 'env-nofilter' if use_env else 'nofilter'
            if use_pairing:
                mode = "pair"
                use_templates = False
                use_env = False

            out_dir = f'{prefix}_{mode}'
            shutil.rmtree(out_dir)
    else:
        a3m_lines = []
        template_paths = []
        for seq in unique_sequences:
            a3m_lines.append(f'>1\n{seq}\n')
            template_paths.append(None)

    # Store MMseqs2 output into dictionary.
    for a3m, templates in zip(a3m_lines, template_paths):
        sequence = a3m.splitlines()[1]
        raw_inputs[sequence] = (a3m, templates)

    # Update with potential custom MSAs.
    for sequence in custom_msas:
        if sequence in raw_inputs:
            templates = raw_inputs[sequence][1]
            raw_inputs[sequence] = (custom_msas[sequence], templates)
        else:
            raw_inputs.update({sequence: (custom_msas[sequence], None)})

    return raw_inputs


def update_msa_lines(old_lines, query_seq, n_pre_gap, n_post_gap):
    new_lines = []
    update = False
    first = True
    for old_line in old_lines:
        if update:
            if first:
                new_lines.append(query_seq)
                update = False
                first = False
            else:
                new_lines.append('-'*n_pre_gap + old_line + '-'*n_post_gap)
                update = False
        else:
            if old_line[0] == '>':
                update = True
            new_lines.append(old_line)
    return new_lines


def runMMseqs2(
        prefix: str,
        sequences: Union[Sequence[str], str],
        use_env: bool = True,
        use_filter: bool = True,
        use_templates: bool = False,
        num_templates: int = 20,
        use_pairing: bool = False,
        host_url: str = 'https://api.colabfold.com'
        ) -> Tuple[Sequence[str], Sequence[Optional[str]]]:
    """ Computes MSAs by querying MMseqs2 API. """

    submission_endpoint = 'ticket/pair' if use_pairing else 'ticket/msa'
    #print(submission_endpoint)
    
    def submit(seqs: Sequence[str], mode: str, N=101) -> Dict[str, str]:
        """ Submits a query of sequences to MMseqs2 API. """

        n, query = N, ''
        for seq in seqs:
            query += f'>{n}\n{seq}\n'
            n += 1

        res = requests.post(f'{host_url}/{submission_endpoint}',
                            data={'q': query, 'mode': mode})
        try:
            out = res.json()
        except ValueError:
            out = {'status': 'ERROR'}

        print(out)
        return out

    def status(ID: int) -> Dict[str, str]:
        """ Obtains the status of a submitted query. """
        res = requests.get(f'{host_url}/ticket/{ID}')
        try:
            out = res.json()
        except ValueError:
            out = {'status': 'ERROR'}

        return out

    def download(ID: int, path: str) -> None:
        """ Downloads the completed MMseqs2 query. """
        res = requests.get(f'{host_url}/result/download/{ID}')
        with open(path, 'wb') as out:
            out.write(res.content)

    # Make input sequence a list if not already.
    sequences = [sequences] if isinstance(sequences, str) else sequences
            
    # Set the mode for MMseqs2.
    if use_filter:
        mode = 'env' if use_env else 'all'
    else:
        mode = 'env-nofilter' if use_env else 'nofilter'
        
    if use_pairing:
        mode = "pair"
        use_templates = False
        use_env = False

    # Set up output path.
    out_path = f'{prefix}_{mode}'
    os.makedirs(out_path, exist_ok=True)
    tar_gz_file = os.path.join(out_path, 'out.tar.gz')
    N, REDO = 101, True

    # Deduplicate and keep track of order.
    unique_seqs = []
    [unique_seqs.append(seq) for seq in sequences if seq not in unique_seqs]
    Ms = [N + unique_seqs.index(seq) for seq in sequences]

    # Call MMseqs2 API.
    if not os.path.isfile(tar_gz_file):
        while REDO:
            # Resubmit job until it goes through
            out = submit(seqs=unique_seqs, mode=mode, N=N)
            while out['status'] in ['UNKNOWN', 'RATELIMIT']:
                # Resubmit
                time.sleep(10 + 5 * (len(unique_seqs) // 50) + random.randint(0, 5))
                out = submit(seqs=unique_seqs, mode=mode, N=N)

            if out['status'] == 'ERROR':
                raise Exception('MMseqs2 API is giving errors. Please confirm '
                                'your input is a valid protein sequence. If '
                                'error persists, please try again in an hour.')

            if out['status'] == 'MAINTENANCE':
                raise Exception('MMseqs2 API is undergoing maintenance. Please '
                                'try again in a few minutes.')
                
            # Wait for job to finish
            ID = out['id']
            while out['status'] in ['UNKNOWN', 'RUNNING', 'PENDING']:
                time.sleep(5 + 5 * (len(unique_seqs) // 50) + random.randint(0, 5))
                out = status(ID)

            if out['status'] == 'COMPLETE':
                REDO = False

            if out['status'] == 'ERROR':
                REDO = False
                raise Exception('MMseqs2 API is giving errors. Please confirm '
                                'your input is a valid protein sequence. If '
                                'error persists, please try again in an hour.')
        # Download results
        download(ID, tar_gz_file)

    # Get and extract a list of .a3m files.
    if use_pairing:
        a3m_files = [os.path.join(out_path, 'pair.a3m')]
    else:
        a3m_files = [os.path.join(out_path, 'uniref.a3m')]
        if use_env:
            a3m_files.append(
                os.path.join(out_path, 'bfd.mgnify30.metaeuk30.smag30.a3m'))
    if not os.path.isfile(a3m_files[0]):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(out_path)

    # Get templates if necessary.
    if use_templates:
        templates = {}
        
        # Read MMseqs2 template outputs and sort templates based on query seq.
        with open(os.path.join(out_path, 'pdb70.m8'), 'r') as f:
            for line in f:
                p = line.rstrip().split()
                M, pdb = p[0], p[1]
                M = int(M)
                if M not in templates:
                    templates[M] = []
                templates[M].append(pdb)

        # Obtain template structures and data files
        template_paths = {}
        for k, TMPL in templates.items():
            TMPL_PATH = os.path.join(prefix+'_'+mode, f'templates_{k}')
            if not os.path.isdir(TMPL_PATH):
                os.mkdir(TMPL_PATH)
                TMPL_LINE = ','.join(TMPL[:num_templates])
                # Obtain the .cif and data files for the templates
                os.system(
                    f'curl -s -L {host_url}/template/{TMPL_LINE} '
                    f'| tar xzf - -C {TMPL_PATH}/')
                # Rename data files
                os.system(
                    f'cp {TMPL_PATH}/pdb70_a3m.ffindex '
                    f'{TMPL_PATH}/pdb70_cs219.ffindex')
                os.system(f'touch {TMPL_PATH}/pdb70_cs219.ffdata')
            template_paths[k] = TMPL_PATH

    # Gather .a3m lines.
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        print(a3m_file)
        with open(a3m_file, 'r') as f:
            for line in f:
                if len(line) > 0:
                    # Replace NULL values
                    if '\x00' in line:
                        line = line.replace('\x00', '')
                        update_M = True
                    if line.startswith('>') and update_M:
                        M = int(line[1:].rstrip())
                        update_M = False
                        if M not in a3m_lines:
                            a3m_lines[M] = []
                    a3m_lines[M].append(line)

    # Return results.
    a3m_lines = [''.join(a3m_lines[n]) for n in Ms]

    if use_templates:
        template_paths_ = []
        for n in Ms:
            if n not in template_paths:
                template_paths_.append(None)
            else:
                template_paths_.append(template_paths[n])
        template_paths = template_paths_
    else:
        template_paths = []
        for n in Ms:
            template_paths.append(None)

    if isinstance(sequences, str):
        return (a3m_lines[0], template_paths[0])
    else:
        return (a3m_lines, template_paths)

def getCustomMSADict(custom_msa_path: str) -> Dict[str, str]:

    custom_msa_dict = {}
    
    onlyfiles = [f for f in os.listdir(custom_msa_path)
                 if os.path.isfile(os.path.join(custom_msa_path, f))]

    custom_msa_dict = {}
    
    for filename in onlyfiles:
        extension = filename.split('.')[-1]
        if extension == 'a3m':
            with open(os.path.join(custom_msa_path, filename)) as f:
                a3m_lines = f.read()
            
            update_seq, seq = True, None
            capture_seq = False
            for line in a3m_lines.splitlines():
                if len(line) > 0:
                    if '\x00' in line:
                        line = line.replace('\x00', '')
                        update_seq = True
                    if line.startswith('>') and update_seq:
                        capture_seq = True
                        update_seq = False
                        header = line
                        continue
                    if capture_seq:
                        seq = line.rstrip()
                        capture_seq = False
                        if seq not in custom_msa_dict:
                            custom_msa_dict[seq] = [header]
                        else:
                            continue

                    if len(line) > 0:
                        custom_msa_dict[seq].append(line)
    
    for seq in custom_msa_dict:
        custom_msa_dict[seq] = '\n'.join(custom_msa_dict[seq])

            #if sequence in custom_msa_dict:
            #    raise ValueError(
            #        f'Multiple custom MSAs found for the sequence the same '
            #        f'sequence: {sequence}. There can only be one custom MSA '
            #        f'per sequence.')

    if custom_msa_dict == {}:
        raise ValueError(
            f'No custom MSAs detected in {custom_msa_path}. Double-check the '
            f'path or no not provide the --custom_msa_path argument. Note that'
            f'custom MSAs must be in .a3m format')
    
    return custom_msa_dict
    

def getMSA(
        sequence: str,
        raw_inputs_from_sequence: RawInput) -> parsers.Msa:

    # Get single-chain MSA.
    raw_inputs = copy.deepcopy(raw_inputs_from_sequence[sequence])
    a3m_lines = raw_inputs[0]
    template_paths = raw_inputs[1]

    return [parsers.parse_a3m(a3m_string=a3m_lines)]


def getUniprotMSA(
        sequence: str,
        raw_inputs_from_sequence: Optional[RawInput] = None,
        ) -> parsers.Msa:
    """ This function essentially creates an MSA with no information. This 
    needs to be updated once Uniprot can be searched with MMseqs2. """

    #logger.warning('AF2 is using an empty UniProt MSA. Results may not be '
    #               'as accurate. This will be changed in the future.')
    
    # Get uniprot MSA
    a3m = f'>{utils.get_hash(sequence)}\n{sequence}\n'

    uniprot_msa = [parsers.parse_a3m(a3m_string=a3m)]

    return uniprot_msa


def getChainFeatures(
        sequences: Sequence[str],
        raw_inputs: RawInput,
        proc_id = None,
        use_templates: bool = False,
        use_multimer = True,
        max_template_date='2100-01-01') -> MutableMapping[str, pipeline.FeatureDict]:
    features_for_chain = {}
    #print(sequences)
    if len(sequences) == 1 or use_multimer:
        for sequence_idx, sequence in enumerate(sequences):
            feature_dict = {}
            # Get sequence features
            feature_dict.update(pipeline.make_sequence_features(
                sequence=sequence, description='query', num_res=len(sequence)))

            # Get MSA features
            msa = getMSA(
                sequence=sequence, raw_inputs_from_sequence=raw_inputs)
            feature_dict.update(pipeline.make_msa_features(msas=msa))

            if len(set(sequences)) > 1:
                uniprot_msa = getUniprotMSA(
                    sequence=sequence)
                valid_feats = msa_pairing.MSA_FEATURES + (
                    'msa_uniprot_accession_identifiers',
                    'msa_species_identifiers',
                )
                all_seq_features = {
                    f'{k}_all_seq': v for
                    k, v in pipeline.make_msa_features(uniprot_msa).items()
                    if k in valid_feats}
                feature_dict.update(all_seq_features)
        
            # Get template features
            if not use_templates:
                feature_dict.update(
                    notebook_utils.empty_placeholder_template_features(
                        num_templates=0, num_res=len(sequence)))

            features_for_chain[
                protein.PDB_CHAIN_IDS[sequence_idx]] = feature_dict
    else:
        feature_dict = {}

        a3m_lines = pair_msa(sequences, raw_inputs)

        total_sequence = ''
        Ls = []
        for sequence in sequences:
            total_sequence += sequence
            Ls.append(len(sequence))

        msa = parsers.parse_a3m(a3m_lines)

        # Sequence features.
        feature_dict.update(
            pipeline.make_sequence_features(
                sequence=total_sequence, description='none', num_res=len(total_sequence)))

        # MSA features.
        feature_dict.update(
            pipeline.make_msa_features([msa]))

        # Template features.
        feature_dict.update(
            notebook_utils.empty_placeholder_template_features(
                num_templates=0, num_res=len(sequence)))

        feature_dict['residue_index'] = chain_break(feature_dict['residue_index'], Ls)
        feature_dict['asym_id'] = np.array(
            [int(n) for n, l in enumerate(Ls) for _ in range(0, l)])

        features_for_chain[
                protein.PDB_CHAIN_IDS[0]] = feature_dict
        
    return features_for_chain


def chain_break(idx_res, Ls, length=200):
    L_prev = 0
    for L_i in Ls[:-1]:
        idx_res[L_prev+L_i:] += length
        L_prev += L_i

    return idx_res


def pair_msa(sequences: Sequence[str], raw_inputs: RawInput) -> parsers.Msa:
    unique_seqs = []
    for seq in sequences:
        if seq not in unique_seqs:
            unique_seqs.append(seq)

    seqs_cardinality = [0]*len(unique_seqs)
    for seq in sequences:
        seq_idx = unique_seqs.index(seq)
        seqs_cardinality[seq_idx] += 1

    unpaired_msas = []
    for seq in unique_seqs:
        unpaired_msas.append(raw_inputs[seq][0])

    return pad_sequences(unpaired_msas, unique_seqs, seqs_cardinality)


def pad_sequences(
        a3m_lines: Sequence[str], query_sequences: Sequence[str],
        query_cardinality: Sequence[int]) -> str:
    _blank_seq = [
        ('-' * len(seq))
        for n, seq in enumerate(query_sequences)
        for _ in range(query_cardinality[n])]

    a3m_lines_combined = []
    pos = 0
    for n, seq in enumerate(query_sequences):
        for j in range(0, query_cardinality[n]):
            lines = a3m_lines[n].split('\n')
            for a3m_line in lines:
                if len(a3m_line) == 0:
                    continue
                if a3m_line.startswith('>'):
                    a3m_lines_combined.append(a3m_line)
                else:
                    a3m_lines_combined.append(
                        ''.join(_blank_seq[:pos] + [a3m_line] + _blank_seq[pos+1:]))
            pos += 1

    return '\n'.join(a3m_lines_combined)
        
            
def getInputFeatures(
        sequences: Sequence[str],
        chain_features: MutableMapping[str, pipeline.FeatureDict],
        min_num_seq: int = 512,
        use_multimer: bool = True,
        ) -> Union[pipeline.FeatureDict,
                   MutableMapping[str, pipeline.FeatureDict]]:

    if len(sequences) == 1 or not use_multimer:
        return chain_features[protein.PDB_CHAIN_IDS[0]]
    else:
        all_chain_features = {}
        for chain_id, features in chain_features.items():
            all_chain_features[
                chain_id] = pipeline_multimer.convert_monomer_features(
                    features, chain_id)

        all_chain_features = pipeline_multimer.add_assembly_features(
            all_chain_features)

        input_features = feature_processing.pair_and_merge(
            all_chain_features=all_chain_features)

        # Pad MSA to avoid zero-size extra MSA.
        return pipeline_multimer.pad_msa(input_features,
                                         min_num_seq=min_num_seq)

