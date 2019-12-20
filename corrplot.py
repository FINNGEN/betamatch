#! /usr/bin/env python3
import pandas as pd,numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from plotnine import *
import argparse,glob,os
import traceback
import os
import re
import math

r = re.compile("x10",re.IGNORECASE)

def main(plot_data, pheno, fields, x_title, y_title, output_name, pval_field=None, p_threshold=None, exp_betas=False):
    if pval_field is not None:
        ## recode unparsable floats
        if plot_data[pval_field].dtype!=np.float64:
            plot_data[pval_field]=plot_data[pval_field].str.replace(" ","").str.replace(u'×',"x").str.replace(r,"e").str.replace("−","-")
        plot_data=plot_data[[fields[0],fields[1], pval_field]].astype("float64")
        #print(plot_data[pval_field])
        plot_data=plot_data.rename(columns={fields[0]:x_title,fields[1]:y_title, pval_field:"pval"})

        if all(np.isnan(plot_data["pval"])):
            ## some don't have p-value defined.... NAs dropped later so set to not filtering 1
            plot_data["pval"] = 1
            p_threshold = None

        if p_threshold is not None:
            print(f"filtering {p_threshold}")
            plot_data = plot_data[ plot_data.pval < p_threshold ]
    else:
        plot_data=data[[fields[0],fields[1]]].astype("float64")
        plot_data=plot_data.rename(columns={fields[0]:x_title,fields[1]:y_title})
    plot_data=plot_data.dropna(axis="index")
    #check for valid data
    #i.e. more than 2 data points
    if plot_data.shape[0]<2:
        print("Too little data!")
        return
    if exp_betas:
        # align to risk allele in 2nd....
        plot_data[[x_title, y_title]] =plot_data[[x_title, y_title]].apply( lambda x: -1*x if x[1]<0 else x, axis=1)
        plot_data[x_title] = np.exp(plot_data[x_title])
        plot_data[y_title] = np.exp(plot_data[y_title])


    [slope,inter,rval,pval,stderr]=stats.linregress( plot_data[x_title], plot_data[y_title])
    x=np.linspace(-100,100,num=plot_data.shape[0])
    y=inter+slope*x
    linedata=pd.DataFrame({x_title:x, y_title:y})
    perf_corr=pd.DataFrame({x_title:x,y_title:x})

    max_val=np.max( np.abs(list(plot_data[x_title])+list(plot_data[y_title]) ) )

    xlim=(-max_val,max_val)
    ylim=(-max_val,max_val)

    if exp_betas:
        xlim=(np.min(plot_data[x_title]),np.max(plot_data[x_title]))
        ylim=(1,np.max(plot_data[y_title]))
        x=np.linspace(xlim[0],xlim[1],num=plot_data.shape[0])
        y=inter+slope*x
        linedata=pd.DataFrame({x_title:x, y_title:y})
        perf_corr=pd.DataFrame({x_title:x,y_title:x})
    #print(slope,inter,rval,pval,stderr)
    #print(corr)
    breaks=[-max_val,0,max_val]
    plot=(ggplot(data=plot_data,mapping=aes(x=x_title,y=y_title))+
        geom_line(data=linedata,mapping=aes(x=x_title,y=y_title),color="#666666" )+
        geom_line(data=perf_corr,mapping=aes(x=x_title,y=y_title),linetype="dashed",color="#888888" )+
        geom_point(color="red",size=0.5)+
        annotate("text",label="R^2: {:.2f}".format(rval**2),x=xlim[0]+(xlim[1]-xlim[0])*0.25,y=ylim[0]+(ylim[1]-ylim[0])*0.85,size=20 )+
        coord_cartesian(xlim=xlim,ylim=ylim)+
        ggtitle(pheno)+
        #scale_x_continuous(breaks=breaks, labels=["{:.2f}".format(-max_val),"0","{:.2f}".format(-max_val)])+
        #scale_y_continuous(breaks=breaks, labels=["{:.2f}".format(-max_val),"0","{:.2f}".format(-max_val)])+
        #scale_x_continuous(limits=xlim)+
        #scale_y_continuous(limits=ylim)+
        theme_minimal()+
        theme(
            axis_ticks_major=element_line(color="black"),
            axis_ticks_minor=None,
            axis_line_x=element_line(color="black"),
            axis_line_y=element_line(color="black")
        )
          )
    #d=plot.draw()
    return plot+theme(figure_size=(6,6))

if __name__=="__main__":
    parser=argparse.ArgumentParser("A utility for plotting correlations from tsv data")
    parser.add_argument("folder",help="data file folder")
    parser.add_argument("--fields",metavar=("field1","field2"),nargs=2,help="column names for the values to be plotted")
    parser.add_argument("--pval_field", default="pval_ext")
    parser.add_argument("--pval_threshold", type=float)
    parser.add_argument("--exp_values", action='store_true')
    parser.add_argument("--x-title",default="x-axis",help="title for x axis")
    parser.add_argument("--y-title",default="y-axis",help="title for y axis")
    parser.add_argument("--out",default="plots.pdf",help="output file name")
    args=parser.parse_args()
    files=glob.glob("{}/*.csv".format(args.folder) )
    plots = []
    for f in files:
        out_fname=os.path.basename(f).split(".")[0] + args.out
        try:
            plot_data=pd.read_csv(f,sep="\t",na_values="-")
            #extract columns
            pheno = os.path.basename(f).split(".")[0]
            print(f"plotting {f}")
            p = main(plot_data, pheno,args.fields,args.x_title,args.y_title,out_fname, args.pval_field, args.pval_threshold, exp_betas=args.exp_values )
            plots.append(p)
        except Exception as e:
            print(f'An exception occurred while plotting {str(e)} \n{traceback.print_exc()}')
    plots=[x for x in plots if x != None]
    save_as_pdf_pages(plots,filename=args.out)
