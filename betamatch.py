#! /usr/bin/python3
import pandas as pd
import tabix
import argparse
import os,subprocess,glob,shlex,re
from subprocess import Popen,PIPE

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


def get_gzip_header(fname):
    """"Returns header for gzipped tsvs, as that is not currently possible using pytabix
    In: file path of gzipped tsv
    Out: header of tsv as a list of column names"""
    gzip_call=shlex.split("gzip -cd {}".format(fname))
    head_call=shlex.split("head -n 1")
    out=[]
    with Popen(gzip_call,stdout=PIPE) as gz:
        with Popen(head_call,stdin=gz.stdout,stdout=PIPE) as hd:
            for line in hd.stdout:
                out.append(line.decode().strip().split("\t"))
    return out[0]

def pytabix(tb,chrom,start,end):
    """Get genomic region from tabixed file
    In: pytabix handle, chromosome, start of region, end of region
    Out: list of variants in region 
    """
    try:
        retval=tb.querys("{}:{}-{}".format(chrom,start,end))
        return retval
    except tabix.TabixError:
        return []

def match_beta(ext_path, fg_summary, info):
    """
    Match beta for external summary variants and our variants
    """
    unified_prefix="unif_"
    
    full_ext_data= pd.read_csv(ext_path,sep="\t",dtype='object')
    full_ext_data=full_ext_data.fillna(value="-")

    mset='^[acgtACGT]+$'
    matchset1=full_ext_data[info[2]].apply(lambda x:bool(re.match(mset,x)))
    matchset2=full_ext_data[info[3]].apply(lambda x:bool(re.match(mset,x)))
    ext_data=full_ext_data[matchset1 & matchset2].copy()
    invalid_ext_data=full_ext_data[~(matchset1 & matchset2)].copy()
    invalid_ext_data["invalid_data"]="YES"

    ext_data[info[2]]=ext_data[info[2]].apply(lambda x:x.upper())
    ext_data[info[3]]=ext_data[info[3]].apply(lambda x:x.upper())
    #load corresponding data from fg file using tabix
    if not os.path.exists("{}.tbi".format(fg_summary)):
        raise FileNotFoundError("Tabix index for file {} not found. Make sure that the file is properly indexed.".format(fg_summary))
    try:
        tabix_handle = tabix.open(fg_summary)
    except tabix.TabixError as e:
        print("An error occurred when opening file {}. Make sure that the file exists and that it is correctly indexed.".format(fg_summary))
        raise
    tmp_lst=[]
    for _,row in ext_data.iterrows():
        tmp_lst=tmp_lst+list(pytabix(tabix_handle,row[info[0]],row[info[1]], row[info[1]] ) )
    header=get_gzip_header(fg_summary)
    summary_data=pd.DataFrame(tmp_lst,columns=header,dtype='object')
    ext_data[info[4]]=pd.to_numeric(ext_data[info[4]])
    summary_data[info[4]]=pd.to_numeric(summary_data[info[4]])
    unif_alt="{}alt".format(unified_prefix)
    unif_ref="{}ref".format(unified_prefix)
    summary_data[[unif_ref,unif_alt]]=summary_data.loc[:,[ info[2], info[3] ]].apply(lambda x: flip_unified_strand(*x),axis=1,result_type="expand")
    ext_data[[unif_ref,unif_alt]]=ext_data.loc[:,[ info[2], info[3] ]].apply(lambda x: flip_unified_strand(*x),axis=1,result_type="expand")
    unif_beta="{}beta".format(unified_prefix)
    summary_data[unif_beta]=summary_data[info[4]]
    ext_data[unif_beta]=ext_data[info[4]]

    for _,row in summary_data.iterrows():
        if list( row[[unif_ref,unif_alt]] )!=sorted( row[[unif_ref,unif_alt]] ):
            summary_data.loc[_,unif_beta]= -1*summary_data.loc[_,unif_beta]
            tmp=summary_data.loc[_,unif_ref]
            summary_data.loc[_,unif_ref]=summary_data.loc[_,unif_alt]
            summary_data.loc[_,unif_alt]=tmp
    for _,row in ext_data.iterrows():
        if list( row[[unif_ref,unif_alt]] )!=sorted( row[[unif_ref,unif_alt]] ):
            ext_data.loc[_,unif_beta]= -1*ext_data.loc[_,unif_beta]
            tmp=ext_data.loc[_,unif_ref]
            ext_data.loc[_,unif_ref]=ext_data.loc[_,unif_alt]
            ext_data.loc[_,unif_alt]=tmp
    ext_data=pd.concat([ext_data,invalid_ext_data],sort=False)
    joined_data=pd.merge(ext_data, summary_data,how="left",on=[info[0],info[1],unif_alt,unif_ref],suffixes=("_ext","_fg"))
    joined_data["beta_same_direction"]=(joined_data["{}_ext".format(unif_beta) ]*joined_data["{}_fg".format(unif_beta) ])>=0
    field_order=["trait", "#chrom", "pos",# "maf", "maf_cases", "maf_controls",  "rsids", "nearest_genes",
     "ref_ext", "alt_ext",  "ref_fg", "alt_fg", 
     "beta_ext", "beta_fg", "se", "sebeta", "pval_ext", "pval_fg",#"unif_ref_ext", "unif_alt_ext", "unif_ref_fg", "unif_alt_fg", 
     "unif_ref","unif_alt","unif_beta_ext", "unif_beta_fg", "invalid_data", "beta_same_direction"]
    joined_data=joined_data[field_order]
    joined_data.to_csv("{}.matched_betas.csv".format(os.path.basename(ext_path).split(".")[0]),index=False,sep="\t",na_rep="-"  )
    return joined_data

def main(args):
    #glob file lists from folders
    ext_files=glob.glob("{}/*.csv".format(args.folder) ) 
    filenames = [os.path.basename(p).split(".")[0] for p in ext_files]
    print(filenames)
    for f in ext_files:
        pheno=os.path.basename(f).split(".")[0]
        print(pheno)
        pheno_path="{}{}.gz".format(args.summaryfolder,pheno)
        out=match_beta(f,pheno_path,args.info)
        #write result to file


if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Match beta of summary statistic and external summaries")
    parser.add_argument("--folder",type=str,help="Folder containing the external summaries that are meant to be used. Files should be names like FinnGen phenotypes.")
    parser.add_argument("--summaryfolder",type=str,help="Finngen summary statistic folder")
    parser.add_argument("--info",nargs=6,default=("#chrom","pos","ref","alt","beta","pval"),metavar=("#chrom","pos","ref","alt","beta","pval"),help="column names")
    args=parser.parse_args()
    main(args)