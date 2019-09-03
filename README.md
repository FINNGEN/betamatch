# betamatch
Tools for matching betas between cherry-picked GWAS studies and FINNGEN results
## betamatch.py
```
usage: betamatch.py [-h] --folder FOLDER --summaryfolder SUMMARYFOLDER --info
                    #chrom pos ref alt beta pval --match-file MATCH_FILE

Match beta of summary statistic and external summaries

optional arguments:
  -h, --help            show this help message and exit
  --folder FOLDER       Folder containing the external summaries that are
                        meant to be used. Files should be names like FinnGen
                        phenotypes.
  --summaryfolder SUMMARYFOLDER
                        Finngen summary statistic folder
  --info #chrom pos ref alt beta pval
                        column names
  --match-file MATCH_FILE
                        List containing the comparisons to be done, as a tsv
                        with columns FG and EXT
```
## corrplot.py
```

usage: A utility for plotting correlations from tsv data [-h]
                                                         [--fields field1 field2]
                                                         [--pval_field PVAL_FIELD]
                                                         [--pval_threshold PVAL_THRESHOLD]
                                                         [--exp_values]
                                                         [--x-title X_TITLE]
                                                         [--y-title Y_TITLE]
                                                         [--out OUT]
                                                         folder

positional arguments:
  folder                data file folder

optional arguments:
  -h, --help            show this help message and exit
  --fields field1 field2
                        column names for the values to be plotted
  --pval_field PVAL_FIELD
  --pval_threshold PVAL_THRESHOLD
  --exp_values
  --x-title X_TITLE     title for x axis
  --y-title Y_TITLE     title for y axis
  --out OUT             output file name
```
