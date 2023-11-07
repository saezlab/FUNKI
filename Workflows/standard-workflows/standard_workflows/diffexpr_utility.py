#from re import S
import sys, yaml, dill, re
#sys.path.insert(1, '/Users/hanna/Documents/projects/Workflows/Python/scUtilities/v4')
#from scUtilities import analysis_params
#import sc_analysis_loops as scl
from copy import deepcopy
from pathlib import Path
from os.path import exists
from os import path, makedirs
import collections
import scanpy as sc, numpy as np, decoupler as dc, matplotlib.pyplot as plt, seaborn as sns, matplotlib as mpl, pandas as pd
from .sc_analysis_baseclass import AnalysisI
from .sc_analysis_baseclass import Baseanalysis


class DiffExpr(AnalysisI):
    """ This class provides methods to run DeSeq2 with bulk data. """
    def __init__(self):
        super().__init__()
        self.paths['diffExpr'] = {}
        self.paths['diffExpr']['deseq2_fig'] = path.join(self.paths['figpath'], 'diffExpr', 'deseq2')
        self.paths['diffExpr']['deseq2_res'] = path.join(self.paths['resultpath'], 'diffExpr', 'deseq2')