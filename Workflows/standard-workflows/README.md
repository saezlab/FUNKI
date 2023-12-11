# Readme

Glossary of shortcuts:
```glossary

ds
    dataset

ap
    analysis_params, config file
    
rawType
    The name for the contant of the raw property of the AnnData object. For example 'counts'.
    
xType
    The name for the contant of the X property of the AnnData object. For example 'filtered'.

ENS
    Ensembl ID

dds
    Eitheer a Deseq2DataSet(AnnData) or a Dds object that holds a Deseq2DataSet as data 

axis
    axis 0 is rows and axis 1 is columns

acts 
    activities (t-values calculated by Decoupler)

    
```


- *preprocessing_utility.py (preu)*
  - Input: bulk dataset
  - plot obs columns
  - get highest expr genes
  - plot violins
  - plot expr filter, filter expr- get highly variable genes
  - plot pcas heatmaps and umaps (for genes and metadata) and metadata associations for various dimensions, get pca loadings
- *diffexpr_utility.py (deu)*
  - Runs DeSeq2 
  - Created a dds class that can hold a DeSeq2DataSet  together with the information about the design factor, reference level and the contrast results.
  - Usage by running ds.run_deseq() per element in ap['diffExpr']['deseq2']['design_factors']. This adds a list of dds objects as attribute 'ddss' to the dataset. Then run for every dds get_contrasts() to save the contrasts in those objects and create the volcano plots.
- *nfcore_utility.py (nfu)*
  -  added option to only rerun necessary parts
- *sc_analysis_baseclass*: metadata paths are added and metadata is read in
- *decoupler_utility (dcu)*
  - model can be saved as csv not just pickle
  - added function for msigdb gmt files to read them in, format the data fitting for decoupler and filter it by substring in collection name (i.e. 'immune')
- *utility_functions*: added otsu threshold function for splitting numerical metadata columns into high/low.
- Features of analysis_params.py: 
  - keyword *take_only*:  
  *take_only* takes a subset of keys after the *project_params* and *dataset_params*.  
  Usage example: *project_params* holds information about all nfCore pipelines that are used in the project. *dataset_params* takes the correct pipeline for every dataset. For example no pipeline for *dataset1* and the *rnaseq* pipeline for *dataset2*:  
    ```
    'dataset1': {   
    'nfcore':{  
        'take_only': []  
    },  
    }, 
    'dataset2': {  
    'nfcore':{
        'pipeline':{
        'take_only': ['rnaseq']  
    }}}
    ```
  
**TestCases**  

  - getpath  
  - merge_dicts  


**Installed packages**   

```glossary

scikit-misc
    ???

PyComplexHeatmap 
    needed for the decoupler functions in get_pca_meta_associations()  
```


**TBD**

- Remove unnecessary packages
- 