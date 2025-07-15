<img style='padding: 10px 10px 20px 10px;' src='assets/funki_logo.svg' width='500'>

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

This package is being developed by Nicolàs Palacio-Escat based on previous work
from Hanna Schumacher and Rosa Hernansaiz-Ballesteros as part of the Saezlab,
Institute for Computational Biomedicine, Univeristy Hospital Heidelberg.

This work is being funded by the CRC 1550 - Molecular Circuits of Heart Disease.

## Preface
Everything within FUNKI is works around an instance of a
[`DataSet`](https://saezlab.github.io/FUNKI/html/input.html#funki.input.DataSet)
class. This object inherits (i.e. is built on top of) the
[`anndata.AnnData`](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html)
class. If you are not familiar with this data structure, you can read more about it their
[publication](https://joss.theoj.org/papers/10.21105/joss.04371) or explore
the package [documentation](https://anndata.readthedocs.io/en/stable/index.html).

## Installation
To install FUNKI in your local computer you can use `pip` as follows:

```bash
pip install git+https://github.com/saezlab/funki
```

You can also install it by cloning or downloading the code in this repository as
described below.

## Usage
If you want to access FUNKI via the graphical interface, you can download the
executable [here](https://github.com/saezlab/FUNKI/releases). If there is no
executable build for your OS, you can also execute the application script with
the Python interpreter as follows:

First you need to download and install FUNKI using the commands below:

```bash
git clone git@github.com:saezlab/FUNKI.git
cd FUNKI
pip install ./
```

Now you can launch the FUNKI interface by simply executing the application
script with:

```bash
python src/funki/app/main.py
```

## Documentation
You can explore all implemented functionalities available
programmatically in the [Documentation](https://saezlab.github.io/FUNKI/).

For GUI usage, you can explore the help pages within the application's "Help"
menu.

---

Developed by Nicolàs Palacio-Escat - Saezlab 2025
