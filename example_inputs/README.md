**This folder holds a couple of possible input files for FUNKI, namely:**   

- CRC1550_Retreat_heartExample.csv  
  -> The dataset that is used for the ChatGPT demonstration.
- Files from the project [GSE50760](https://github.com/adugourd/GSE50760/tree/main/example_inputs):  
  - Differential_genes.csv  
    -> Only a list of genes that uses ora. There is no difference between using a file with a genes column or entering the genes in the textfield on the website. 
  - Differential_pvals.csv  
    -> Genes and pvals column, runs with ulm 
  - Differential_stats.csv  
    -> Genes and t vals column, runs with ulm, is used for the testdata option in FUNKI
      
      
**Please make sure that input csv files are formatted like this:** 

- Comma separated
- Numbers use the english notation with dots instead of commas
- First column are genes
- The optional second column holds statistics (LFC, pval, ...) 
- Columns have names
- There can be more columns but they are ignored by FUNKI
