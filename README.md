<img style='padding: 10px 10px 20px 10px;' src='./docs/source/_images/funki_logo.svg' width='500'>

# Welcome to FUNKI
Welcome to FUNKI, the omics FUNctional analysis worKflows Interface tool. This
Python package is intended to integrate different omic data analysis workflows
including a graphical user interface (GUI), but also as a standalone Python
package that users can integrate into their existing pipelines.

## Disclaimer
This package is currently in development and features are being added regularly.
If you have any ideas/suggestions or if you find any bug, please feel free to
open an [GitHub issue](https://github.com/saezlab/FUNKI/issues). You can also
contribute via [pull request](https://github.com/saezlab/FUNKI/pulls).

Please note that the GUI is not yet implemented and FUNKI can only be accessible
as a Python package.

## Installation
To install FUNKI in your local computer you can use `pip` as follows:

```bash
pip install git+https://github.com/saezlab/funki
```

## Usage
Everything within FUNKI is works around an instance of a
[`DataSet`](https://saezlab.github.io/FUNKI/html/input.html#funki.input.DataSet)
object. This class in turn inherits (i.e. is built on top of) the
[`anndata.AnnData`](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html)
class

## Documentation
You can find the full documentation in the
[FUNKI documentation webpage](https://saezlab.github.io/FUNKI/).
