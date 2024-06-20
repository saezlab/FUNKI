Welcome to FUNKI, the omics FUNctional analysis worKflows Interface tool. This
Python package is intended to integrate different omic data analysis workflows
including a graphical user interface (GUI), but also as a standalone Python
package that users can integrate into their existing pipelines.

You are currently using the GUI, if you are interested in the Python
package, please visit the [FUNKI GitHub page](https://github.com/saezlab/FUNKI)

### Usage

You can navigate the different sections through the tabs on the left panel.
Further information on the different inputs and parameters can be obtained by
hovering your mouse over the red question marks accompanying each of them.

- **Data:** Here you will be able to upload your data and metadata of your
  samples as well as visualize the distribution of values and categories. 
  Currently supported formats include `.csv`, `.txt` and `.xlsx` but more are
  compatible formats are planned in the future. The data provided is assumed to
  contain the observations (e.g. cells, samples, etc.) on the rows and the 
  variables (e.g. genes, proteins, annotations, etc.) on the columns. Currently
  we assume genes are defined by Gene Symbols (currently only needed for the
  enrichment step), but we are planning to add ID conversion tools in the near
  future.
- **Filter & normalization:** This part of FUNKI is currently designed mainly
  for single-cell transcriptomic data sets, but support for other types of omics
  data are being worked on. In this section you will be able to select among
  several filters for minimum and maximum genes per cell as well as the maximum
  percentage of mitochondrial genes per cell (assumes the gene symbols are 
  preceded by the label `'MT-'`). Here you can also normalize your counts based 
  on a target number of counts per cell and also log-transform the data. Several
  quality control plots are shown below your selection and are updated as soon
  as the changes are applied to the data.
- **Clustering:** On this section you can apply unsupervised clustering by using
  community detection algorithms. Currently implemented are the two most common
  approaches: Louvain and Leiden (which is a modified version of Louvain aiming
  to improve Louvain's issues with disconnected communities). In this panel you
  can also visualize your data through common dimensionality reduction 
  embeddings, including Principal Component Analysis (PCA), t-Stochastic 
  Neighbor Embedding (tSNE) and Unifold Manifold Approximation and Projection
  (UMAP). You can then choose to color your observations based on the results of
  the clustering methods or any other metadata variable provided.
- **Enrichment:** Here you will be able to select a gene set collection from
  a comprehensive list of resources providing molecular information. The list 
  includes major databases, resources and publications like MSigDB, KEGG, 
  CellPhoneDB, PROGENy among others. Once the resource is selected, you can 
  further specify your serach by applying filters based on the different 
  information contained in the given resource. Once you have the gene set 
  collection of interest, enrichment analysis can then be applied.