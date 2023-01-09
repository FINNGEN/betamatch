#! /usr/bin/env python3
import pandas as pd,numpy as np #type: ignore
import plotnine as p9
from beta_utils import *


def plot_simple(datapoints:pd.DataFrame,regression_data:RegressionResults,weighted_regression_data:RegressionResults,x_title:str,y_title:str,r2:float,weighted_r2:float,title:str)->p9.ggplot:
    """Plot the correlation plot
    Args:
        datapoints: point data (DataFrame) with columns x,y, x_se, y_se
        regression_data: regression data for non-weighted regression
        weighted_regression_data: regression data for weighted regression
        x_title: x title
        y_title: y title
        r2: float
        weighted_r2: float
        title: plot title
    """
    x_ci = 1.96*datapoints["x_se"]
    y_ci = 1.96*datapoints["y_se"]
    datapoints["ci_x_neg"] = datapoints["x"]-x_ci
    datapoints["ci_x_pos"] = datapoints["x"]+x_ci
    datapoints["ci_y_neg"] = datapoints["y"]-y_ci
    datapoints["ci_y_pos"] = datapoints["y"]+y_ci

    min_val = min(np.min(datapoints["x"].values),np.min(datapoints["y"].values))
    max_val = max(np.max(datapoints["x"].values),np.max(datapoints["y"].values))
    xlim=(-np.abs(min_val),max_val)
    ylim=(-np.abs(min_val),max_val)
    plot = (
        p9.ggplot(data=datapoints,mapping=p9.aes(x="x",y="y"))+
        p9.geom_point(color="red",size=0.5)+
        p9.geom_errorbar(mapping=p9.aes(x="x",y="y",ymin = "ci_y_neg",ymax="ci_y_pos"))+
        p9.geom_errorbarh(mapping=p9.aes(x="x",y="y",xmin = "ci_x_neg",xmax="ci_x_pos"))+
        p9.geom_abline(intercept=0,slope=1,color="#666666")+
        p9.geom_abline(intercept=0,slope=regression_data.slope,color="#888888",linetype="dashed")+
        p9.geom_abline(intercept=0,slope=weighted_regression_data.slope,color="#666666",linetype="dashdot")+
        p9.annotate("text",label="R^2 (pearson):{:>5.2g}  slope:{:>5.2g}  se(slope):{:>5.2g}".format(r2,regression_data.slope,regression_data.slope),x=xlim[0]+(xlim[1]-xlim[0])*0.5,y=ylim[0]+(ylim[1]-ylim[0])*0.99,size=10 )+
        p9.annotate("text",label="R^2 (weighted):{:>5.2g}  slope:{:>5.2g}  se(slope):{:>5.2g}".format(weighted_r2,weighted_regression_data.slope,weighted_regression_data.stderr),x=xlim[0]+(xlim[1]-xlim[0])*0.5,y=ylim[0]+(ylim[1]-ylim[0])*0.94,size=10 )+
        p9.coord_cartesian(xlim=xlim,ylim=ylim)+
        p9.ggtitle(title)+
        p9.xlab(x_title)+
        p9.ylab(y_title)+
        p9.theme(
            axis_ticks_major=p9.element_line(color="black"),
            axis_ticks_minor=None,
            axis_line_x=p9.element_line(color="black"),
            axis_line_y=p9.element_line(color="black"),
            figure_size=(6,6)
        )
    )
    return plot