# Super-fast-Manhattan-plots
Super fast plotting of Manhattan of GWAS. Don't spend your time waiting but plotting! Thanks to Yanick's magic R code


RUN: `Rscript generic_manhattan.R GWASdata.tsv chrom pos minus_log10_p TRUE outputfile 0.00000005 snps.to.color.tsv 1000000 black FALSE`

- 1st = filename
- 2,3,4th = column names of chrom, pos and P 
- 5th = TRUE = if Pval is logtransformed or not. 
- 6th=outputfile
- 7th = P Treshold
- 8th = a file (e.g. 'snps.to.color.tsv' includes position of top-snp to color, use following if you don't want to color anything:
```
      | chromosome | position | colour |
      |------------|----------|--------|
      | 19         | 0        | grey   |
```
- 9th = window of loci to be coloroued (e.g. 1000000 bases)
- 10th = color of ? not sure anymore  ()
- 11th = custom ylim, if FALSE then take the maximum, otherwise put a number in like 20, so max = -logp=20 


Code is not cleaned in tidy funcions, appologies. 
