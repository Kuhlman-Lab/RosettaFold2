""" Utility functions for parsing queries from a list of input files. """

# Standard imports.
import os
import csv
import logging
import random
from typing import Sequence, Tuple, Union

# Constants.
restypes = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V'
]

# Get logger.
logger = logging.getLogger('query_utils')

# (filename, sequence)
MonomerQuery = Tuple[str, str]

# (filename, oligomer_state, [sequences])
MultimerQuery = Tuple[str, str, Sequence[str]]

# (filename, [sequences])
CleanQuery = Tuple[str, Sequence[str]]

# List of 20 standard amino acids.
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def generate_random_sequences(
        lengths: Sequence[Union[int, Sequence[int]]],
        num_seq: int,
        aalist: Sequence[str] = None) -> Sequence[Sequence[str]]:

    if aalist == None:
        aalist = AA_LIST

    seqs_list = []
    for length in lengths:
        for _ in range(num_seq):
            seqs = []
            if isinstance(length, int):
                seqs.append(''.join(random.choices(aalist, k=length)))
            elif isinstance(length, list):
                for chain_length in length:
                    seq = ''.join(random.choices(aalist, k=chain_length))
                    seqs.append(seq)
            seqs_list.append(seqs)
        
    return seqs_list

def parse_fasta(fasta_string: str) -> Tuple[Sequence[str], Sequence[str]]:
  """Parses FASTA string and returns list of strings with amino-acid sequences.

  Arguments:
    fasta_string: The string contents of a FASTA file.

  Returns:
    A tuple of two lists:
    * A list of sequences.
    * A list of sequence descriptions taken from the comment lines. In the
      same order as the sequences.
  """
  sequences = []
  descriptions = []
  index = -1
  for line in fasta_string.splitlines():
    line = line.strip()
    if line.startswith('>'):
      index += 1
      descriptions.append(line[1:])  # Remove the '>' at the beginning.
      sequences.append('')
      continue
    elif not line:
      continue  # Skip blank lines.
    sequences[index] += line

  return sequences, descriptions

def parse_fasta_files(files: Sequence[str]) -> Sequence[MonomerQuery]:
    """ 
    Parse a list of .fasta files and return a list of monomer queries. 
    """
    query_list = []
    
    for filename in files:
        with open(filename, 'r') as f:
            fasta_string = f.read()

        seqs, _ = parse_fasta(fasta_string)

        # .fasta files can contain multiple query sequences.
        for seq in seqs:
            query_list.append( (filename, seq) )

    return query_list


def parse_csv_files(files: Sequence[str]) -> Sequence[MultimerQuery]:
    """ 
    Parse a list of .csv files and return a list of multimer queries. 
    """
    query_list = []

    for filename in files:
        with open(filename, newline='') as f:
            reader = csv.reader(f, delimiter=',')
            # Each row is a single query with possibly many sequences
            # summarized by the oligomer state.
            for row in reader:
                oligomer = row[0]
                og_sequences = row[1:]
                sequences = []
                for sequence in og_sequences:
                    if '#' in sequence:
                        sequences.append(sequence.split('#')[0])
                    else:
                        if sequence:
                            sequences.append(sequence)

                query_list.append( (filename, oligomer, sequences) )

    return query_list


def clean_and_validate_queries(
        input_queries: Union[Sequence[MonomerQuery], Sequence[MultimerQuery]],
        min_length: int,
        max_length: int,
        max_multimer_length: int
    ) -> Sequence[CleanQuery]:
    """ 
    Validates and cleans input queries. 
    """
    query_list = []

    for query in input_queries:
        query = _clean_and_validate_single_query(
            query=query,
            min_length=min_length,
            max_length=max_length,
            max_multimer_length=max_multimer_length)
        query_list.append(query)
        
    if len(query_list) > 0:
        return query_list
    
    else:
        raise ValueError('No files contain a valid query, please provide at '
                         'least one query.')

    
def _clean_and_validate_single_query(
        query: Union[MonomerQuery, MultimerQuery], min_length: int,
        max_length: int, max_multimer_length: int) -> CleanQuery:
    """
    Checks that the parsed query is ok and returns a clean version of it.
    """
    filename = query[0]
    sequences = query[-1]
    
    if filename == '_INPUT_':
        oligomer = ''
    else:
        if len(query) == 2:
            # If a monomer query is given, then it has an oligomeric state of 1.
            oligomer = '1'
        else:
            oligomer = query[1]

    # If a monomer query is given, then it has single sequence. Need to treat as
    # a list of sequences.
    if isinstance(sequences, str):
        sequences = [sequences]
    
    # Clean filename by removing all parent directories.
    clean_filename = os.path.basename(filename)

    # Remove whitespaces, tabs, and end lines and uppercase all sequences.
    clean_sequences = []
    for sequence in sequences:
        clean_sequence = sequence.translate(
            str.maketrans('', '', ' \n\t')).upper()
        aatypes = set(restypes) # 20 canonical aatypes.
        if not set(clean_sequence).issubset(aatypes):
            raise ValueError(
                f'Query parsed from {clean_filename} has a sequence with '
                f'non-amino acid letters: {set(clean_sequence) - aatypes}. '
                f'AlphaFold only supports 20 standard amino acids as inputs.')
        if len(clean_sequence) < min_length:
            raise ValueError(
                f'Query parsed from {clean_filename} has a sequence that is '
                f'too short: {len(clean_sequence)} amino acids, while the '
                f'minimum is {min_length}.')
        if len(clean_sequence) > max_length:
            raise ValueError(
                f'Query parsed from {clean_filename} has a sequence that is '
                f'too long: {len(clean_sequence)} amino acids, while the '
                f'maximum is {max_length}. If you believe you have the '
                f'resources for this query, overwrite the default max_length '
                f'by providing the argument: --max_length NEW_MAX.')
        clean_sequences.append(clean_sequence)

    if len(clean_sequences) < 1:
        raise ValueError(
            f'Query parsed from {clean_filename} does not have any detectable '
            f'sequences.')

    # Clean oligomer and validate shape
    if oligomer == '':
        if filename != '_INPUT_':
            logger.warning(f'Inferring oligomeric state from sequences provided in '
                           f'{clean_filename}.')
        clean_oligomer = ':'.join(['1'] * len(clean_sequences))
    else:
        clean_oligomer = oligomer.translate(
            str.maketrans('', '', ' \n\t'))

    oligomer_vals = set('123456789:')
    if not set(clean_oligomer).issubset(oligomer_vals):
        raise ValueError(
            f'Query parsed from {clean_filename} has an oligomer state '
            f'with non-valid characters: '
            f'{set(clean_oligomer) - oligomer_vals}.')
    oligos = clean_oligomer.split(':')
    if len(oligos) > len(clean_sequences):
        raise ValueError(
            f'Query parsed from {clean_filename} has more oligomeric '
            f'states than number of sequences: {len(oligos)} > '
            f'{len(clean_sequences)}. Oligomer is {clean_oligomer}.')
    if len(oligos) < len(clean_sequences):
        raise ValueError(
            f'Query parsed from {clean_filename} has less oligomeric '
            f'states than number of sequences: {len(oligos)} < '
            f'{len(clean_sequences)}. Oligomer is {clean_oligomer}.')

    clean_sequences = [[seq] * int(oligo)
        for seq, oligo in zip(clean_sequences, oligos)]
    clean_sequences = [chain for homomer in clean_sequences for chain in homomer]
    
    total_multimer_length = sum([len(chain) for chain in clean_sequences])
    if total_multimer_length > max_multimer_length:
        raise ValueError(
            f'Query parsed from {clean_filename} has a total multimer length '
            f'that is too long: {total_multimer_length}, while the maximum '
            f'is {max_multimer_length}. If you believe you have the resources '
            f'for this query, overwrite the default max_multimer_length by '
            f'providing the argument: --max_multimer_length NEW_MAX.')
    elif total_multimer_length > 1536:
        logger.warning(f'The accuracy of the multimer system has not been '
                       f'fully validated above 1536 residues. Query from '
                       f'{clean_filename} is a total length of '
                       f'{total_multimer_length}.')
    
    # If there is only one sequence, then the query is a monomer query.
    # If there is more than one sequence, then the query is a multimer query.
    return (clean_filename, clean_sequences)

    
def detect_duplicate_queries(
        query_list: Sequence[CleanQuery]) -> Sequence[CleanQuery]:
    """ 
    Detects duplicate queries from query list. If a same query comes from 
    two different sources, it is considered a duplicate. 
    """
    clean_query_list = []

    for query in query_list:
        # If the clean_query_list is empty, query is not a dupe.
        if len(clean_query_list) == 0:
            clean_query_list.append(query)
        else:
            dupe = False

            # If the sequences of the queries are the same, then its a dupe.
            for old_query in clean_query_list:
                if old_query[1] == query[1]:
                    dupe = True

            if dupe == False:
                clean_query_list.append(query)

    # Sort the clean_query_list such that monomer queries appear first.
    clean_query_list = sorted(clean_query_list, key=lambda x: len(x[1]))

    return clean_query_list
    

def getFullSequence(query: Union[MonomerQuery, MultimerQuery]) -> str:
    """
    Given a query, returns the full sequence for which a structure is predicted.
    I.e. if given a multimer, returns the full multimer sequence.
    """
    if len(query) == 2:
        full_sequence = query[1]
    else:
        oligomer = query[1]
        sequences = query[2]

        oligo_list = oligomer.split(':')

        full_sequence = ''.join([
            seq * int(oligo) for seq, oligo in zip(sequences, oligo_list)])

    return full_sequence
