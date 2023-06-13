**This folder holds a couple of possible input files for FUNKI, namely:**   

- CRC1550_Retreat_heartExample.csv  
  -> The dataset that is used for the ChatGPT demonstration.
- differential_genes.csv  
  -> Only a list of genes that uses ora. There is no difference between using a file with a genes column or entering the genes in the textfield on the website. 
- differential_pvals.csv  
  -> genes and pvals column, runs with ulm
- differential_stats.csv  
  -> genes and t vals column, runs with ulm  
      
      
**Please make sure that input csv files are formatted like this:** 

- comma separated
- numbers use the english notation with dots instead of commas
- first column are genes
- the optional second column holds statistics (LFC, pval, ...) 
- columns have names
- there can be more columns but they are ignored by FUNKI
