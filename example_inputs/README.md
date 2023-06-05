example_inputs/differential_stats.csv:
  the stats can be fold changes, t-values or zscores. They range between - and + infinity.
  run decoupleR's ulm
example_inputs/differential_pvals.csv
  the stats are p-values. They range between 0 and 1
  use all the genes to generate a background similar to background_genes.csv (see below)
  use a p-value threshold of < 0.05 and generate an input identical to the differential_genes.csv file (see below)
  run decoupleR's ORA using the differential genes and back gournd genes as input
example_inputs/differential_genes.csv
  run decoupleR's ORA using the default example_inputs/background_genes.csv file as background genes
example_inputs/background_genes.csv
  background genes to use as default background of ORA
  
CHECK WITH PAU IF DECOUPLER ORA CAN BE RUN WITH INPUTS list OF SIGNIFICANT GENES AND BACKGORUND



Where do the tables come from?
Genelist -> field for # of background genes