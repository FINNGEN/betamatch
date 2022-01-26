#! /usr/bin/python3
from typing import NamedTuple
import pandas as pd, numpy as np
import tabix
import argparse
import os,subprocess,glob,shlex,re
from subprocess import Popen,PIPE
from scipy.stats import pearsonr, norm

from beta_utils import *

class FGCols(NamedTuple):
    chr:str
    pos:str
    ref:str
    alt:str
    beta:str
    pval:str
    se:str

class ExtCols(NamedTuple):
    chr:str
    pos:str
    ref:str
    alt:str
    beta:str
    pval:str
    se:str
    study_doi:str

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

def flip_beta(r1,a1,beta1,beta2):
    """If beta1 <0, flip betas (and consequently alleles)
    """
    if beta1 < 0:
        return (a1,r1,-beta1,-beta2)
    else:
        return (r1,a1,beta1,beta2)

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
        return list(retval)
    except tabix.TabixError:
        return []


def match_beta(ext_path, fg_summary, info_ext:ExtCols, info_fg:FGCols):
    """
    Match beta for external summary variants and our variants
    In: ext fpath, fg fpath, column tuple
    Out: df containing the results. DOES NOT SAVE FILES
    """
    unified_prefix="unif_"
    ext_dtype = {info_ext.chr:object,
                info_ext.pos:object,
                info_ext.ref:object,
                info_ext.alt:object,
                info_ext.beta:float,
                info_ext.pval:float,                
                info_ext.se:float,
                info_ext.study_doi:object}
    fg_dtype = {info_fg.chr:object,
                info_fg.pos:object,
                info_fg.ref:object,
                info_fg.alt:object,
                info_fg.beta:float,
                info_fg.pval:float,
                info_fg.se:float}
    full_ext_data= pd.read_csv(ext_path,sep="\t",dtype=ext_dtype)
    full_ext_data[[info_ext.ref,info_ext.alt]]=full_ext_data[[info_ext.ref,info_ext.alt]].fillna(value="-")

    full_ext_data[info_ext.pval]=full_ext_data[info_ext.pval].astype(float)
    full_ext_data[info_ext.beta]=full_ext_data[info_ext.beta].astype(float)
    #replace missing se values with values derived from beta+pvalue
    for idx,row in full_ext_data.iterrows():
        if pd.isna(row[info_ext.se]):
            #calculate se
            zscore=np.abs(norm.ppf(row[info_ext.pval]/2) )
            se=np.abs(row[info_ext.beta])/zscore
            full_ext_data.loc[idx,info_ext.se] = np.nan if se <= 0 else se
    full_ext_data[info_ext.se]=full_ext_data[info_ext.se].astype(float)
    mset='^[acgtACGT]+$'
    matchset1=full_ext_data[info_ext.ref].apply(lambda x:bool(re.match(mset,x)))
    matchset2=full_ext_data[info_ext.alt].apply(lambda x:bool(re.match(mset,x)))
    ext_data=full_ext_data[matchset1 & matchset2].copy()
    invalid_ext_data=full_ext_data[~(matchset1 & matchset2)].copy()
    invalid_ext_data["invalid_data"]="YES"
    ext_data[info_ext.ref]=ext_data[info_ext.ref].apply(lambda x:x.upper())
    ext_data[info_ext.alt]=ext_data[info_ext.alt].apply(lambda x:x.upper())
   
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
        tmp_lst=tmp_lst+pytabix(tabix_handle,row[info_ext.chr],row[info_ext.pos],row[info_ext.pos] )
    header=get_gzip_header(fg_summary)

    summary_data=pd.DataFrame(tmp_lst,columns=header).astype(dtype=fg_dtype)
    ext_data[info_ext.beta]=pd.to_numeric(ext_data[info_ext.beta],errors='coerce')
    summary_data[info_fg.beta]=pd.to_numeric(summary_data[info_fg.beta])
    unif_alt="{}alt".format(unified_prefix)
    unif_ref="{}ref".format(unified_prefix)
    summary_data[[unif_ref,unif_alt]]=summary_data.loc[:,[ info_fg.ref, info_fg.alt ]].apply(lambda x: flip_unified_strand(*x),axis=1,result_type="expand")
    ext_data[[unif_ref,unif_alt]]=ext_data.loc[:,[ info_ext.ref, info_ext.alt ]].apply(lambda x: flip_unified_strand(*x),axis=1,result_type="expand")
    unif_beta="{}beta".format(unified_prefix)
    summary_data[unif_beta]=summary_data[info_fg.beta]
    ext_data[unif_beta]=ext_data[info_ext.beta]

    summary_data["sort_dir"] = summary_data[[unif_ref,unif_alt]].apply(lambda x: -1 if (sorted(list(x) ) != list(x)) else 1,axis=1)
    summary_data[[unif_ref,unif_alt]]=summary_data[[unif_ref,unif_alt]].apply(lambda x: sorted(list(x)),axis=1,result_type="expand")
    summary_data[unif_beta]=summary_data["sort_dir"]*summary_data[unif_beta]
    summary_data=summary_data.drop(labels="sort_dir",axis="columns")
    ext_data["sort_dir"] = ext_data[[unif_ref,unif_alt]].apply(lambda x: -1 if (sorted(list(x)) != list(x)) else 1,axis=1)
    ext_data[[unif_ref,unif_alt]]=ext_data[[unif_ref,unif_alt]].apply(lambda x: sorted(list(x)),axis=1,result_type="expand")
    ext_data[unif_beta]=ext_data["sort_dir"]*ext_data[unif_beta]
    ext_data=ext_data.drop(labels="sort_dir",axis="columns")
    ext_data=pd.concat([ext_data,invalid_ext_data],sort=False)
    info_fg_rename = {info_fg[i]:info_ext[i] for i in range(len(info_ext)-1)} 
    summary_data.rename(columns=info_fg_rename, inplace=True)
    joined_data=pd.merge(ext_data, summary_data,how="left", on=[info_ext[0],info_ext[1],unif_alt,unif_ref],suffixes=("_ext","_fg"))

    unif_beta_ext="{}_ext".format(unif_beta)
    unif_beta_fg="{}_fg".format(unif_beta)
    joined_data[[unif_ref,unif_alt,unif_beta_ext,unif_beta_fg]]=joined_data[[unif_ref,unif_alt,unif_beta_ext,unif_beta_fg]].apply(
        lambda x: flip_beta(*x),axis=1,result_type="expand")
    
    joined_data["beta_same_direction"]=(joined_data["{}_ext".format(unif_beta) ]*joined_data["{}_fg".format(unif_beta) ])>=0
    if "trait" not in joined_data.columns:
        traitcol = []
    else:
        traitcol = ["trait"]
    field_order=traitcol + [info_ext.chr, info_ext.pos,
     info_ext.ref+"_ext",info_ext.alt+"_ext", 
     info_ext.ref+"_fg", info_ext.alt+"_fg", 
     info_ext.beta+"_ext",info_ext.beta+"_fg",
     info_ext.pval+"_ext",info_ext.pval+"_fg",
     info_ext.se+"_ext",info_ext.se+"_fg",
     "unif_ref","unif_alt","unif_beta_ext", "unif_beta_fg", "invalid_data", "beta_same_direction", info_ext.study_doi]

    joined_data=joined_data[field_order]
    return joined_data

def extract_doi(joined_data, info):
    """
    Output doi strings
    In: joined data, column in ext file with doi
    Out: string with all dois contained in the file
    """
    doi_concat=','.join(joined_data[info].dropna().unique())
    return doi_concat

def main(info_ext:ExtCols,info_fg:FGCols,match_file,out_f,pval_filter):
    """
    Match betas between external summ stats and FG summ stats
    In: folder containing ext summaries, folder containing fg summaries, column tuple, matching tsv file path 
    Out:  
    """
    match_df=pd.read_csv(match_file,sep="\t",header=None,names=["EXT","FG"])
    output_list=[]
    r2s=[]
    for _,row in match_df.iterrows():
        ext_path = row["EXT"]
        fg_path = row["FG"]
        fg_name = os.path.splitext(os.path.basename(fg_path))[0]
        ext_name = os.path.splitext(os.path.basename(ext_path))[0]
        #check existance
        print(row)
        output_fname="{}x{}.betas.tsv".format(ext_name.split(".")[0],fg_name)
        if (os.path.exists( ext_path ) ) and ( os.path.exists( fg_path ) ):
            matched_betas=match_beta(ext_path,fg_path,info_ext,info_fg)
            matched_betas=matched_betas[matched_betas[info_ext.pval+"_ext"]<=pval_filter]
            stat_data=matched_betas[["unif_beta_ext","unif_beta_fg",info_ext.se+"_ext"]].dropna(axis="index",how="any")
            dois_ext=extract_doi(matched_betas,info_ext.study_doi)
            if not stat_data.empty:
                r2,w_r2,n_r,n_w=calculate_r2(stat_data,"unif_beta_ext","unif_beta_fg",info_ext.se+"_ext")
                normal_regression = calculate_regression(stat_data["unif_beta_ext"].values,stat_data["unif_beta_fg"].values )
                weighted_regression = calculate_regression(stat_data["unif_beta_ext"].values,stat_data["unif_beta_fg"].values,1/(stat_data[info_ext.se+"_ext"]**2) )
                row={"phenotype":output_fname.split(".")[0],"R^2":r2,"Weighted R^2 (1/ext var)":w_r2,"N (unweighted)":n_r,"N (weighted)":n_w, "study_doi": dois_ext}
                row.update( {"Regression slope":normal_regression.slope,"Weighted regression slope":weighted_regression.slope,"Regression intercept":0.0,
                    "Weighted regression intercept":0.0,
                    "Regression std.err.":normal_regression.stderr,
                    "Weighted regression std.err.":weighted_regression.stderr,
                    "Regression slope p-value": normal_regression.pval,
                    "Weighted regression slope p-value": weighted_regression.pval} )
                r2s.append(row)

            matched_betas.to_csv(path_or_buf=out_f+"/"+output_fname,index=False,sep="\t",na_rep="-")
            output_list.append(output_fname)
        else:
            print("One of the files {}, {} does not exist. That pairing is skipped.".format(ext_path,fg_path))
    r2s=pd.DataFrame(r2s)
    r2s.to_csv("r2_table.tsv",sep="\t",index=False,float_format="%.3g",na_rep="-")
    print("The following files were created:")
    [print(s) for s in output_list]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Match beta of summary statistic and external summaries")
    #parser.add_argument("--folder",type=str,required=True,help="Folder containing the external summaries that are meant to be used. Files should be names like FinnGen phenotypes.")
    #parser.add_argument("--summaryfolder",type=str,required=True,help="Finngen summary statistic folder")
    parser.add_argument("--info-ext",nargs=8,required=True,default=("#chrom","pos","ref","alt","beta","pval","se","study_doi"),metavar=("#chrom","pos","ref","alt","beta","pval","se","study_doi"),help="column names for external file")
    parser.add_argument("--info-fg",nargs=7,required=True,default=("#chrom","pos","ref","alt","beta","pval","se"),metavar=("#chrom","pos","ref","alt","beta","pval","se"),help="column names for finngen file")
    parser.add_argument("--match-file",required=True,help="List containing the comparisons to be done, as a tsv with columns FG and EXT")
    parser.add_argument("--output-folder",required=True,help="Output folder")
    parser.add_argument("--pval-filter",default=1.0,type=float,help="Filter p-value for summary file")
    args=parser.parse_args()
    extcols = ExtCols(
        args.info_ext[0],
        args.info_ext[1],
        args.info_ext[2],
        args.info_ext[3],
        args.info_ext[4],
        args.info_ext[5],
        args.info_ext[6],
        args.info_ext[7]
    )
    fgcols = FGCols(
        args.info_fg[0],
        args.info_fg[1],
        args.info_fg[2],
        args.info_fg[3],
        args.info_fg[4],
        args.info_fg[5],
        args.info_fg[6]
    )

    main(extcols,fgcols,args.match_file,args.output_folder,args.pval_filter)
