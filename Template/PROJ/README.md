# PROJ  

Collaboration Form: *insert link to collaboration form*

**How to use different parameter settings**  
  
The analysis parameters are defined in *analysis_params.py*. After changing the file, run it with *poetry run* or execute it from within VSCode.  


## Analysis Steps   

- Bulk 
  - nf-core rnaseq pipeline
  - preprocessing
    - highly expressed genes
    - DeSeq2
      - Decoupler
        - progeny
        - collecri
        - MSigDB
      - Liana
  - more filtering
    - highly variable genes (hvg)
    - pca (loadings, loadings heatmap, umap per gene and metadata column)
- Single Cell
  - Decoupler (acts + umaps, mean acts + heatmap)

## Installations

**PyEnsembl**  
For translating the Ensembl IDs into gene symbols  

> !pip install pyensembl
> from pyensembl import EnsemblRelease
> # Install the needed release:  
> !pyensembl install --release 109 --species mouse

### Citations  
If this pipeline was used, the following tools should be cited: 

- Decoupler, Omnipath, CollecTri, Progeny, MSigDB
- Scanpy/AnnData
- nf-core rnaseq pipeline
- Standard-Workflows repo