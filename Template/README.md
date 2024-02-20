# Template Usage

You can use the 'Proj'-folder as template for a new project. 

The template includes the basic folder structure for a bulk dataset which we call in this case '01', a single cell dataset called '01' and a pseudobulk dataset called '01_sc'. You can choose different names. 

**Needed changes/preparations**  

*In general:*  

- Clone the FUNKI repository to get the template and the standard-workflows package
- Rename 'PROJ' to the shortcut of your project
- Install the poetry environment in the scripts folder
- Adjust the analysis_params script and execute it

*Bulk:*  
- Add your fastq files in the data folder
- Adjust the metadata tables
- Download the MSigDB gmt files that you want to include in the analysis. Save them in the priorKnowledge folder.
- Run the main_bulk script
  - If you start with raw data, you can use the script to prepare the run of the rnaseq nf-core pipeline. Then run the pipeline. Afterwards, start again from the beginning of the script.

*SingleCell:*  
- Add single cell count tables per sample or a h5ad file to the data folder
- Run the main_sc script step by step. You might need to make some adjustments. 
- You now have a Pseudobulk data file. Move it to the 01_sc data folder and run the main_01_sc script.