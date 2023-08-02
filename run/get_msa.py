import os, sys, argparse
import numpy as np
from collections import OrderedDict, Counter
from string import ascii_uppercase, ascii_lowercase
import hashlib, re, os

sys.path.append('/proj/kuhl_lab/RosettaFold2/RoseTTAFold2/')
from input_prep.api import run_mmseqs2
from input_prep.setup import QueryManager

alphabet_list = list(ascii_uppercase+ascii_lowercase)

def get_hash(x): 
  return hashlib.sha1(x.encode()).hexdigest()

def get_unique_sequences(seq_list):
    unique_seqs = list(OrderedDict.fromkeys(seq_list))
    return unique_seqs
  
def getCustomMSADict(custom_msa_path: str):

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

def calc_msa(seq, jobname, cov=50, id=90, max_msa=2048,
            mode="unpaired_paired"):

  assert mode in ["unpaired","paired","unpaired_paired"]
  seqs = [seq] if isinstance(seq,str) else seq

  # collapse homooligomeric sequences
  counts = Counter(seqs)
  u_seqs = list(counts.keys())
  u_nums = list(counts.values())

  # expand homooligomeric sequences
  first_seq = "/".join(sum([[x]*n for x,n in zip(u_seqs,u_nums)],[]))
  msa = [first_seq]

  path = os.path.join(jobname,"msa")
  os.makedirs(path, exist_ok=True)
  if mode in ["paired","unpaired_paired"] and len(u_seqs) > 1:
    print("getting paired MSA")
    out_paired = run_mmseqs2(u_seqs, f"{path}/", use_pairing=True)
    headers, sequences = [],[]
    for a3m_lines in out_paired:
      n = -1
      for line in a3m_lines.split("\n"):
        if len(line) > 0:
          if line.startswith(">"):
            n += 1
            if len(headers) < (n + 1):
              headers.append([])
              sequences.append([])
            headers[n].append(line)
          else:
            sequences[n].append(line)
    
    # filter MSA
    with open(f"{path}/paired_in.a3m","w") as handle:
      for n,sequence in enumerate(sequences):
        handle.write(f">n{n}\n{''.join(sequence)}\n")
    os.system(f"hhfilter -i {path}/paired_in.a3m -id {id} -cov {cov} -o {path}/paired_out.a3m")
    with open(f"{path}/paired_out.a3m","r") as handle:
      for line in handle:
        if line.startswith(">"):
          n = int(line[2:])
          xs = sequences[n]
          # expand homooligomeric sequences
          xs = ['/'.join([x]*num) for x,num in zip(xs,u_nums)]
          msa.append('/'.join(xs))
  
  if len(msa) < max_msa and (mode in ["unpaired","unpaired_paired"] or len(u_seqs) == 1):
    print("getting unpaired MSA")
    out = run_mmseqs2(u_seqs,f"{path}/")
    Ls = [len(seq) for seq in u_seqs]
    sub_idx = []
    sub_msa = []
    sub_msa_num = 0
    for n,a3m_lines in enumerate(out):
      sub_msa.append([])
      with open(f"{path}/in_{n}.a3m","w") as handle:
        handle.write(a3m_lines)
      # filter
      os.system(f"hhfilter -i {path}/in_{n}.a3m -id {id} -cov {cov} -o {path}/out_{n}.a3m")
      with open(f"{path}/out_{n}.a3m","r") as handle:
        for line in handle:
          if not line.startswith(">"):
            xs = ['-'*l for l in Ls]
            xs[n] = line.rstrip()
            # expand homooligomeric sequences
            xs = ['/'.join([x]*num) for x,num in zip(xs,u_nums)]
            sub_msa[-1].append('/'.join(xs))
            sub_msa_num += 1
      sub_idx.append(list(range(len(sub_msa[-1]))))
    
    while len(msa) < max_msa and sub_msa_num > 0:
      for n in range(len(sub_idx)):
        if len(sub_idx[n]) > 0:
          msa.append(sub_msa[n][sub_idx[n].pop(0)])
          sub_msa_num -= 1
        if len(msa) == max_msa:
          break

  with open(f"{jobname}/msa.a3m","w") as handle:
    for n,sequence in enumerate(msa):
      handle.write(f">n{n}\n{sequence}\n")
      
def get_msa(sequence, 
            jobname = "test", #@param {type:"string"}
            order = 1, #@param ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"] {type:"raw"}
            msa_method = "mmseqs2", #@param ["mmseqs2","single_sequence","custom_a3m"]
            pair_mode = "unpaired_paired", #@param ["unpaired_paired","paired","unpaired"] {type:"string"}
            sym = "X", #@param ["X","C", "D", "T", "I", "O"]
            collapse_identical = False,#@param {type:"boolean"}
            max_msa = 256, #@param [16, 32, 64, 128, 256, 512] {type:"raw"}
            custom_msa_path=None):
  
  # process
  max_extra_msa = max_msa * 8

  sequence = re.sub("[^A-Z:]", "", sequence.replace("/",":").upper())
  sequence = re.sub(":+",":",sequence)
  sequence = re.sub("^[:]+","",sequence)
  sequence = re.sub("[:]+$","",sequence)

  if sym in ["X","C"]:
    copies = order
  elif sym in ["D"]:
    copies = order * 2
  else:
    copies = {"T":12,"O":24,"I":60}[sym]
    order = ""
  symm = sym + str(order)

  sequences = sequence.replace(":","/").split("/")
  if collapse_identical:
    u_sequences = get_unique_sequences(sequences)
  else:
    u_sequences = sequences
  sequences = sum([u_sequences] * copies,[])
  lengths = [len(s) for s in sequences]

  sequence = "/".join(sequences)
  #jobname = jobname+"_"+symm+"_"+get_hash(sequence)[:5]

  print(f"jobname: {jobname}")
  print(f"lengths: {lengths}")

  os.makedirs(jobname, exist_ok=True)
  if msa_method == "mmseqs2":
    calc_msa(u_sequences, jobname, mode=pair_mode, max_msa=max_extra_msa)

  elif msa_method == "single_sequence":
    u_sequence = "/".join(u_sequences)
    with open(f"{jobname}/msa.a3m","w") as a3m:
      a3m.write(f">{jobname}\n{u_sequence}\n")

  #TODO: custom_msa testing
  elif msa_method == "custom_a3m":
    if not custom_msa_path:
      raise ValueError("Please provide a path to a custom a3m file")
    with open(custom_msa_path,"r") as handle:
      lines = handle.readlines()
    a3m_lines = []
    for line in lines:
      line = line.replace("\x00","")
      if len(line) > 0 and not line.startswith('#'):
        a3m_lines.append(line)

    with open(f"{jobname}/msa.a3m","w") as a3m:
      a3m.write("\n".join(a3m_lines))
      
def main(args):
  #parse sequences
  qm = QueryManager(
        input_dir=args.input_dir,
        #sequences=[args.sequence],
        sequences=[],
        min_length=args.min_length,
        max_length=args.max_length,
        max_multimer_length=args.max_multimer_length)
  qm.parse_files()
  qm.parse_sequences()

  queries = qm.queries
  print(f'Queries have been parsed. {len(queries)} queries found.')

  #for each sequence, get msa
  for query, i in zip(queries, range(len(queries))):
    sequence = ":".join(query[1])
    get_msa(sequence, 
            jobname = os.path.join(args.output_dir, "query_"+str(i)),
            msa_method=args.msa_mode,
            pair_mode=args.pair_mode,
            collapse_identical=args.collapse_identical,
            sym=args.symmetry,
            order=args.order,
            custom_msa_path=args.custom_msa_path)


if __name__ == "__main__":
    
    # Parse arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--input_dir', 
                        type=str, 
                        default='inputs/',
                        help='Input directory where .fasta or .csv files live.')

    parser.add_argument('--msa_mode', type=str, 
                        default='mmseqs2', 
                        choices=["mmseqs2", "single_sequence", "custom_a3m"], 
                        help='Mode to generate MSAs.')
    
    parser.add_argument('--custom_msa_path', type=str, 
                        default=None, 
                        help='If using custom_a3m mode, path to custom a3m file.')
    
    parser.add_argument('--output_dir', 
                        type=str, 
                        default='mmseqs2/', 
                        help='Where to write the output .a3m MSA files.')
    
    parser.add_argument('--pair_mode', 
                        default='unpaired_paired', 
                        choices=["unpaired_paired","paired","unpaired"], 
                        help='unpaired - generate seperate MSA for each protein.'
                        'unpaired_paired - attempt to pair sequences from the same genome.'
                        'paired - only use sequences that were successfully paired.')
    
    parser.add_argument('--symmetry', 
                        default='X', 
                        choices=["X","C", "D", "T", "I", "O"], 
                        help='[C]yclic, [D]ihedral, [T]etrahedral, [I]cosahedral, [O]ctahedral and [X] for unknown.')
    
    parser.add_argument('--order', 
                        default=1, 
                        choices=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"], 
                        help='define number of copies in X/C/D symmetries')
    
    parser.add_argument('--collapse_identical', 
                        action='store_true',
                        help='Default is False. Remove identical (homooligomeric) sequences from input sequence. '
                        '(This could be useful if you plan to use the symmetry option)')
    
    #not implemented yet
    parser.add_argument('--min_length', type=int, default=50, help='Minimum length of MSA matches returned.')
    parser.add_argument('--max_length', type=int, default=10000, help='Maximum length of MSA matches returned.')
    parser.add_argument('--max_multimer_length', type=int, default=25000, help='Maximum length of MSA multimer matches.')

    args = parser.parse_args()
    
    # Get MSAs
    main(args)