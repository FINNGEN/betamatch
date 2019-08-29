#! /usr/bin/python3
import pandas as pd,numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from plotnine import *
import argparse,glob,os

def main(fname, fields, x_title, y_title, output_name):
    #load data
    data=pd.read_csv(fname,sep="\t",na_values="-")
    #extract columns
    plot_data=data[[fields[0],fields[1]]].astype("float64")
    plot_data=plot_data.rename(columns={fields[0]:x_title,fields[1]:y_title})
    plot_data=plot_data.dropna(axis="index")
    [slope,inter,rval,pval,stderr]=stats.linregress(np.array(plot_data).T)
    x=np.linspace(-100,100,num=plot_data.shape[0])
    y=inter+slope*x
    linedata=pd.DataFrame({x_title:x, y_title:y})
    perf_corr=pd.DataFrame({x_title:x,y_title:x})
    max_val=np.max( np.abs(list(plot_data[x_title])+list(plot_data[y_title]) ) )
    xlim=(-max_val,max_val)
    ylim=(-max_val,max_val)
    #xlim=(np.min(plot_data[x_title]),np.max(plot_data[x_title]) )
    #ylim=(np.min(plot_data[y_title]),np.max(plot_data[y_title]) )
    print(slope,inter,rval,pval,stderr)
    #print(corr)
    breaks=[-max_val,0,max_val]
    plot=(ggplot(data=plot_data,mapping=aes(x=x_title,y=y_title))+
        geom_line(data=linedata,mapping=aes(x=x_title,y=y_title),color="#666666" )+
        geom_line(data=perf_corr,mapping=aes(x=x_title,y=y_title),linetype="dashed",color="#888888" )+
        geom_point(color="red",size=0.5)+
        annotate("text",label="R^2: {:.2f}".format(rval**2),x=xlim[0]+(xlim[1]-xlim[0])*0.25,y=ylim[0]+(ylim[1]-ylim[0])*0.85,size=20 )+
        coord_cartesian(xlim=xlim,ylim=ylim)+
        scale_x_continuous(breaks=breaks, labels=["{:.2f}".format(-max_val),"0","{:.2f}".format(-max_val)])+
        scale_y_continuous(breaks=breaks, labels=["{:.2f}".format(-max_val),"0","{:.2f}".format(-max_val)])+
        theme_minimal()+
        theme(
            axis_ticks_major=element_line(color="black"),
            axis_ticks_minor=None,
            axis_line_x=element_line(color="black"),
            axis_line_y=element_line(color="black")
        )
          )
    #d=plot.draw()
    save_as_pdf_pages([plot+theme(figure_size=(6,6))],filename=output_name)

if __name__=="__main__":
    parser=argparse.ArgumentParser("A utility for plotting correlations from tsv data")
    parser.add_argument("folder",help="data file folder")
    parser.add_argument("--fields",metavar=("field1","field2"),nargs=2,help="column names for the values to be plotted")
    parser.add_argument("--x-title",default="x-axis",help="title for x axis")
    parser.add_argument("--y-title",default="y-axis",help="title for y axis")
    parser.add_argument("--out",default="output",help="output file name")
    args=parser.parse_args()
    files=glob.glob("{}/*.csv".format(args.folder) )
    for f in files:
        out_fname=os.path.basename(f).split(".")[0] + args.out
        try:
            main(f,args.fields,args.x_title,args.y_title,out_fname)
        except:
            pass