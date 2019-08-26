#! /usr/bin/python3
import pandas as pd
import argparse
import os,subprocess,glob

def flip_unified_strand(a1,a2):
    """
    Flips alleles to the A strand if necessary 
    """
    allele_dict={"T":"A","C":"G","G":"C"}
     # check if the A variant is present
    if 'A' not in a1 + a2 and 'a' not in a1+a2:
        # for both/ref and alt map them to the A strand and order each one lexicographically
        a1 = ''.join([allele_dict[elem.upper()] for elem in a1])
        a2 = ''.join([allele_dict[elem.upper()] for elem in a2])
    return (a1,a2)


def match_beta(ext_path, fg_summary, info):
    """
    Match beta for external summary variants and our variants
    """
    #load files

    #drop those that do not have sufficient information, e.g. alleles

    #match alleles, flip betas for those that need it

    #output file

    return None

def main(args):
    #glob file lists from folders
    ext_files=glob.glob("{}/*.csv".format(args.folder) ) 
    filenames = [os.path.basename(p).split(".")[0] for p in ext_files]
    for f in ext_files:
        pheno=os.path.basename(f).split(".")[0]
        pheno_path="{}/{}.gz".format(args.summaryfolder,pheno)
        out=match_beta(f,pheno_path,args.info)
        #write result to file


if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Match beta of summary statistic and external summaries")
    parser.add_argument("--folder",type=str,help="Folder containing the external summaries that are meant to be used. Files should be names like FinnGen phenotypes.")
    parser.add_argument("--summaryfolder",type=str,help="Finngen summary statistic folder")
    parser.add_argument("--info",nargs=6,default=("chr","pos","ref","alt","beta","pval"),metavar=("chr","pos","ref","alt","beta","pval"),help="column names")
    args=parser.parse_args()
    main(args)