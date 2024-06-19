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
class. You can explore all implemented functionalities available
programmatically in the [Documentation](https://saezlab.github.io/FUNKI/).

If you want to access FUNKI via the GUI, you need to execute the application 
script with your Python interpreter, but first you need to download and install
FUNKI using the commands below:

```bash
git clone git@github.com:saezlab/FUNKI.git
cd FUNKI
pip install ./
```

Now you can launch the FUNKI by simply running the application script:

```bash
python src/funki/app.py
```

This should automatically open the application in your default internet browser.
If that is not the case, you can type or copy the following address (default)
in the  browser bar: `http://127.0.0.1:8050/`

**Note:** Refreshing the page will also restart the application.

## Documentation
You can find the full documentation in the
[FUNKI documentation webpage](https://saezlab.github.io/FUNKI/).
