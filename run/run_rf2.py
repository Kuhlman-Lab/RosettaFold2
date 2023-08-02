import torch
import os, sys

sys.path.append('/proj/kuhl_lab/RosettaFold2/RoseTTAFold2/')
from network.predict import namedtuple, Predictor, read_index, read_data

def get_args():
    default_model = "/proj/kuhl_lab/RosettaFold2/RoseTTAFold2/network/weights/RF2_apr23.pt"

    import argparse
    parser = argparse.ArgumentParser(description="RoseTTAFold2NA")
    
    #Recommended to start with this option.
    parser.add_argument("-input_dir", 
                        help="To loop over multiple queries, use this option.", default=None)
    
    #Not recommended unless you know what you're doing.
    parser.add_argument("-inputs", help="R|Input data in format A:B:C, with\n"
                        "   A = multiple sequence alignment file\n"
                        "   B = hhpred hhr file\n"
                        "   C = hhpred atab file\n"
                        "Spaces seperate multiple inputs.  The last two arguments may be omitted\n"
                        "WARNING: TEMPLATES/CUSTOM MSA NOT SET UP YET.\n",
                        default=None, 
                        nargs='+')
    
    #Defaults to None since templates have not been set up yet.
    parser.add_argument("-db", 
                        help="HHpred database location", 
                        default=None)
    
    parser.add_argument("-prefix", 
                        default="S", 
                        type=str, 
                        help="Output file prefix [S]")
    
    parser.add_argument("-symm", 
                        default="C1", 
                        help="Symmetry group (Cn,Dn,T,O, or I).  If provided, MSA should cover the asymmetric unit. [C1]")
    
    parser.add_argument('-msa_concat_mode', 
                        default="diag", 
                        choices=["diag", "repeat", "default"], 
                        help='Default is diag. Defines how the msa is concatenated.'
                        '| repeat | diag | default |'
                        '| AAA    | A--  | A--     |'
                        '|        | -A-  | -AA     |'
                        '|        | --A  |         |')
    
    #right now, the only option is the default weights.
    parser.add_argument("-model", 
                        default=default_model, 
                        help="Model weights. [weights/RF2_apr23.pt]")
    
    parser.add_argument("-n_recycles", 
                        default=3, 
                        type=int, 
                        help="Number of recycles to use [3].")
    
    parser.add_argument("-n_models", 
                        default=1, 
                        type=int, 
                        help="Number of models to predict [1].")
    
    parser.add_argument("-subcrop", 
                        default=-1, 
                        type=int, 
                        help="Subcrop pair-to-pair updates. For very large models (>3000 residues) a subcrop of 800-1600 can improve structure accuracy and reduce runtime. A value of -1 means no subcropping. [-1]")
    
    parser.add_argument("-nseqs", 
                        default=256, 
                        type=int, 
                        help="The number of MSA sequences to sample in the main 1D track [256].")
    
    parser.add_argument("-nseqs_full", 
                        default=2048, 
                        type=int, 
                        help="The number of MSA sequences to sample in the wide-MSA 1D track [2048].")
    
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()

    if (args.db is not None):
        FFDB = args.db
        FFindexDB = namedtuple("FFindexDB", "index, data")
        ffdb = FFindexDB(read_index(FFDB+'_pdb.ffindex'),
                         read_data(FFDB+'_pdb.ffdata'))
    else:
        ffdb = None
        
    print ("Using weights at ", args.model)

    if (torch.cuda.is_available()):
        print ("Running on GPU")
        pred = Predictor(args.model, torch.device("cuda:0"))
    else:
        print ("Running on CPU")
        pred = Predictor(args.model, torch.device("cpu"))
        
    outdir = 'rf2out'
    os.makedirs(outdir, exist_ok=True)
    
    if args.input_dir:
        input_dir = args.input_dir
        queries = os.listdir(input_dir)
        print(queries)
        
        inputs = [os.path.join(input_dir, q) for q in queries if os.path.isdir(os.path.join(input_dir, q))]
        print(inputs)
        input_msa_files = [os.path.join(q, "msa.a3m") for q in inputs if os.path.exists(os.path.join(q, "msa.a3m"))]
        print(input_msa_files)

        for file in input_msa_files:
            outdir_query = os.path.join('rf2out', os.path.basename(os.path.dirname(file)))
            os.makedirs(outdir_query, exist_ok=True)
            pred.predict(
                inputs=[file], 
                out_prefix=outdir_query + "/" + args.prefix,
                symm=args.symm, 
                n_recycles=args.n_recycles, 
                n_models=args.n_models, 
                subcrop=args.subcrop, 
                nseqs=args.nseqs, 
                nseqs_full=args.nseqs_full, 
                ffdb=ffdb)
    
    else:
        pred.predict(
            inputs=args.inputs, 
            out_prefix=outdir + "/" + args.prefix,
            symm=args.symm, 
            n_recycles=args.n_recycles, 
            n_models=args.n_models, 
            subcrop=args.subcrop, 
            nseqs=args.nseqs, 
            nseqs_full=args.nseqs_full, 
            ffdb=ffdb)