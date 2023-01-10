# betamatch
Tools for matching betas between cherry-picked GWAS studies and FINNGEN results
## betamatch.py
```
usage: betamatch.py [-h] --info-ext #chrom pos ref alt beta pval se study_doi --info-fg #chrom pos ref alt beta pval se --match-file MATCH_FILE --output-folder
                    OUTPUT_FOLDER [--pval-filter PVAL_FILTER] [--drop-extra-cols]

Match beta of summary statistic and external summaries

optional arguments:
  -h, --help            show this help message and exit
  --info-ext #chrom pos ref alt beta pval se study_doi
                        column names for external file
  --info-fg #chrom pos ref alt beta pval se
                        column names for finngen file
  --match-file MATCH_FILE
                        List containing the comparisons to be done, as a tsv with columns FG and EXT
  --output-folder OUTPUT_FOLDER
                        Output folder
  --pval-filter PVAL_FILTER
                        Filter p-value for summary file
  --drop-extra-cols     Drop extra columns
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
