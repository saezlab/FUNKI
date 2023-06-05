"""Docstring for the scfunctions.py module.

v1: matches to decoupler_v1, putting the functions into classes

Modules names should have short, all-lowercase names.  The module name may
have underscores if this improves readability.

"""
# for markdown magic
from __future__ import absolute_import
from IPython.core.getipython import get_ipython
from IPython.core.magic import (Magics, magics_class,  cell_magic)
import sys
from io import StringIO
from markdown import markdown
from IPython.core.display import HTML

from pathlib import Path
#from scUtilities import scfunctions_v1 as scfuncs
import sys, itertools, logging.config, os, logging #, h5py, igraph, pandas, louvain, random, re, dill, functools, collections
import scanpy as sc, matplotlib.pyplot as plt, seaborn as sns, numpy as np

#import pytest


from functools import reduce
from copy import deepcopy
from mergedeep import merge

def print_paths(paths):
    """Prints the folder structure that the paths define together with the labels.
    The bigree.print_tree method can only handle paths that have the same root node. 
    Therefore the paths are grouped by their root and then each group is printed. 
    Paths that appear more than once get reduced to one by joining their values together. 

    Args:
        paths (dict): A dictionary with paths as values and the name of the path as key. 
    Example: 
        > d = {'Users/hanna/docs/hu': 'hu', 'Users/hanna/docs/hu/miau': 'huhumiau', 'Users/hanna/docs/hui': 'hiui'}
    Utility: 
        dict.fromkeys(paths, {"": 90} # creates a dictionary with the same keys and the given value.
    """
    import bigtree
    import pandas
    #hm = copy.deepcopy(paths)
    #paths = {k: paths[k] for k in paths} # switch keys and values as the path must be the key
    # paths.loc[paths["key"].startswith("/"), "key"] 

    paths = pandas.DataFrame({'key': paths.values(), '': paths.keys()})
    paths['key'] = '/' + paths['key'].replace('^/', '', regex=True).replace('/$', '', regex=True)
    paths = paths.groupby(['key'])[''].apply(', '.join).reset_index()
    paths['keyStart'] = paths['key'].str.split('/').str[1]
    grp=paths.groupby('keyStart')

    temp_out = StringIO()
    sys.stdout = temp_out


    for k, item in grp:
        p = bigtree.dataframe_to_tree(grp.get_group(k))
        bigtree.print_tree(p, attr_list=[""], attr_bracket=["[", "]"])

    sys.stdout = sys.__stdout__
    return temp_out.getvalue()

def deep_get(dictionary, *keys):
    return reduce(lambda d, key: d.get(key) if d else None, *keys, dictionary)

def dict_delete_key(dict, path):
    """Delete a key from a nested dictionary.

    Args:
        dict (dict): The dictionary where the key shall be removed. The dict is changed in place.
        path (tuple): The path to the key that shall be removed (use getpath() to get the path)
    """    
    pathlen = len(path)
    r = range(0,pathlen-1)
    del (reduce(lambda d,i: d[path[i]], r, dict))[path[pathlen-1]]


def getpath(nested_dict, value, prepath=(), search_value = True) -> tuple:
    """ Get a tuple of keys that lead to a given value in a given nested dictionary.

    Args:
        nested_dict (dict): A nested dictionary that shall be searched through.
        value (char): Either a key (search_value = False) or a value in the dict. 
        prepath (tuple, optional): Internal parameter used for recursion. Defaults to ().
        search_value (bool, optional): If the given 'value' to search for is a key or a value. Defaults to True.

    Returns:
        path (tuple): The order of keys that leads to the value or key in the dict (when searching for a key the key is inclusive in the result path).
    
    Use Case:
        workflows.scfunctions.getpath({1: {'name': 'John', 'age': '27', 'sex': 'Male'}, 2: {'name': 'Marie', 'age': '22', 'sex': 'Female'}}, 'Female')

    Caution: 
        When the value is a key, another key is the same if it is exactly the same or if its prefix (before underscore) is the same. 
    
    Code mostly taken from:
        [Stackoverflow - answer from jsf](https://stackoverflow.com/questions/22162321/search-for-a-value-in-a-nested-dictionary-python)
    """
    
    for k, v in nested_dict.items():
        path = prepath + (k,)
        if(search_value == False): 
            # search for key instead of value
            if (k == value) | (k.split('_')[0] == value): # to allow for using multiple similar keys in a dict like 'replacement_1', ...
                return path
            elif hasattr(v, 'items'):
                p = getpath(v, value, path, search_value) # recursive call
                if p is not None:
                    return p
        elif value in v: # found value
            return path
        elif hasattr(v, 'items'): # v is a dict
            p = getpath(v, value, path, search_value) # recursive call
            if p is not None:
                return p

def dict_replace(dict, v, path) -> dict:
    """ Fill in a value *v* into a dictionary *dict* at the position given by *path* """
    pathlen = len(path)
    r = range(0,pathlen-1)
    (reduce(lambda d,i: d[path[i]], r, dict))[path[pathlen-1]] = v
    return dict

def dict_delete_key(dict, path) -> dict:
    pathlen = len(path)
    r = range(0,pathlen-1)
    del (reduce(lambda d,i: d[path[i]], r, dict))[path[pathlen-1]]

def merge_dicts(dict_1: dict, dict_2: dict) -> dict:
    """
    First updates dict_2 with dict_1. 
    Then, inserts the contents of dict_1 at the positions of the 'add' entries in the values of dict_2. 
    Restrictions: Can handle onyl one level beneath and then only when dict_2 is one level beneath dict_1.

    UseCase
    -------
    dict_1 with default values for an analysis. dict_2 with intended changes to these default values. 
    These can be replacements or additions. Values that should be extended must be placed in a list in both dicts. 
    dict_2 has an 'add' entry wherever the default values from dict_1 should be inserted. 
        
    Parameters
    ----------
    dict_1 : dict
    dict_2 : dict

    Variables & Functions
    ---------------------
    placeholder : 'add'
    keys: all keys till one key value contains the placeholder
    value: the list that contains the placeholder
    ind: index of the placeholder in the list
    getpath() 
    deep_get()
    copy:deepcopy()

    long_var_name : {'hi', 'ho'}, optional
        Choices in brackets, default first when optional.

    Returns
    -------
    merged : dict

    Test Thoughts
    -------------
    Test Cases
    - key in dict_1 but not in dict_2 -> dog isplayful
    - key in dict_2 but not in dict_1 -> cat sound
    - 'add' entry at the beginning, in the middle, at the end -> beginning & middle
    - none, one and multiple 'add' entries -> 2
    - 'add' outside of list
    - 'add' in nested list or tuple value or in not nested dict

    Pseudocode:
        1. Find position of placeholder in dict_2. -> list of keys that lead there
        2. If dict_1 first level keys don't include dict_2 or next iter dict_2:
             a. Get path in dict_1 for key 'next iter dict_2'
             b. Remove last element from keypath
             c. Get value from keypath
             2.1 If placeholder not in dict_2:
                   a. 
    Examples
    --------
    dict_1 can start one layer above dict_2. 
    >>> dict_1 = {'pets': {'dog': {'name': 'Bello', 'sound': ['wuff', 'wau'], 'isplayful': True}, 'cat': {'name': 'Kitty', 'food': ['fish']}}}
    >>> dict_2 = {'dog': {'name': 'Wauwau', 'sound': ['wau', 'add', 'grrr']}, 'cat': {'name': 'Miezy', 'sound': 'maunz', 'food': ['add', 'mice']}}
    >>> merge_dicts(dict_1, dict_2)
    {'pets': {'dog': {'name': 'Wauwau', 'sound': ['wau', 'wuff', 'grrr'], 'isplayful': True}, 'cat': {'name': 'Miezy', 'food': ['fish', 'mice'], 'sound': 'maunz'}}}
    """
    placeholder = 'ADD'
    keys = getpath(dict_2, placeholder)
    #merged = deepcopy(dict_1)
    #merged.update(dict_2)
    if(dict_2 and (next(iter(dict_2)) not in dict_1.keys())):
        p = getpath(deepcopy(dict_1), next(iter(dict_2)), search_value=False)
        p = p[:len(p)-1]
        d = deep_get(dict_1,p)
        if(keys == None):
            merge(d, merge(deepcopy(d), deepcopy(dict_2)))
            return deepcopy(dict_1)
    else: 
        d = dict_1
    if(keys != None): 
        value = list(deep_get(dict_2, keys))
        ind = value.index(placeholder)
        value = list(set(value[:ind]+ deep_get(d, keys) + value[ind+1:]))
        merged = dict_replace(deepcopy(dict_2), value, keys)
        return merge_dicts(d, merged)
    else: # no keys and no level differences
        return merge(deepcopy(dict_1), deepcopy(dict_2)) # to get the missing keys from dict_1    



def replace_dictvalues(target, placeholder = 'REPLACEMENTS') -> dict:
    """ Executes replacements indicated by the placeholder key (or a key matching'placeholder_.*') in a dictionary. 

    Args:
        target (dict): A dictionary. 
        placeholder (str): If the placeholder key is not inside the dictionary, nothing happens.  
        Otherwise the key value pairs in the placeholder dictionary replace all values of the keys that follow.
        Replacements are done in dictionaries on the same level as the placeholder that follow in the ordered dict after the placeholder
        The dictionaries are only checked for having the key in the top level. 
    Returns:
        dict: A copy of the input dictionary with executed replacements
    Usage: 
        Can be used together with (before or after) merge_dicts()
    Example: 
    >>> dict_2 = {'pets': {'snail': {'name': 'Schnucki',  'food': ['Löwenzahn', 'water'], 'sound': ['knusper',]},
    >>>                    'replacements_1': {'name': 'LittleWauWau','sound': ['kikerikii']},    
    >>>                    'dog': {'sound': ['wau', 'add', 'grrr']}, 
    >>>                    'replacements_2': {'color': 'orange'},
    >>>                    'cat': {'name': 'Miezy', 'food': ['add', 'mice'], 'color': 'pink & yellow stripes', 'child': {
    >>>                         'babyCat': {'name': 'MiniMiez', 'food': ['milk'], 'color': 'yellow with pink dots'}}}}}
    >>> replace_dictvalues(dict_2)
    {'pets': {'snail': {'name': 'Schnucki', 'food': ['Löwenzahn', 'water'], 'sound': ['knusper']},
                'dog': {'sound': ['kikerikii'], 'name': 'LittleWauWau'},
                'cat': {'name': 'LittleWauWau', 'food': ['add', 'mice'], 'color': 'orange', , 'sound': ['kikerikii'], 'child': {
                    'babyCat': {'name': 'MiniMiez', 'food': ['milk'], 'color': 'yellow with pink dots'}}}}}
    """
    replacements_path = getpath(target, placeholder, search_value = False)          # ['pets', 'replacements']
    replacements_path_fordeletion = replacements_path
    if replacements_path != None: 
        placeholder = replacements_path_fordeletion[-1]
        replacements_path = list(replacements_path)
        replacement_values = deep_get(target, replacements_path)            # {'name': 'LittleWauWau', 'sound': 'kikerikii'}
        del replacements_path[-1]                                                    # ['pets']  
        replacements_dict = deep_get(target, replacements_path)             # {'replacements':{}, 'dog':{}, 'cat':{}}
        keys = list(replacements_dict.keys())
        ind = keys.index(placeholder)                                                # 0
        for k in keys[ind+1:]:
            if k != placeholder:
                sub_dict_path = getpath(target, k, search_value = False)    # ('pets', 'dog')
                sub_dict = deep_get(target, sub_dict_path)
                merged = merge_dicts(sub_dict, replacement_values)
                target = dict_replace(deepcopy(target), merged, replacements_path+[k]) # cp
        dict_delete_key(target, replacements_path_fordeletion)
        return replace_dictvalues(target)
    else: 
        return target


def add_nested_key(dict, path:list):
    """Adds a key to a dictionary even when the path to that key doesn't already exist

    Args:
        dict (_type_): _description_
        path (list): full path to the key
    """
    if len(path)>0:
        key = path[0]
        path.pop(0)
        if key not in dict:
            add = {key: {}}
            dict.update(add)
        add_nested_key(dict[key], path)


def deleteVariables(varStart) -> None:
   def get_variables(n) -> list :

      return [k for k in globals() if k.startswith(n)]

   vars_for_del = get_variables(varStart)
   for v in vars_for_del:
      del globals()[v]


# look for example at seurat_clusters, should louvain/leiden be added?
def plot_scUmap (datasets) : 
    def plot (dataset, name) : 
        sc.settings.figdir = paths["figpath"] / name 
        sc.pl.umap(dataset, color = dataset.obs, save = '_c12', show = False, components='all') #legend_loc="on data")

        return [plot(datasets[key], key) for key in datasets]

# from Tarrasch "https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists"
def dict_product (dicts):
    """
    >>> list(dict_product(dict(number=[1,2], character='ab')))
    [{'character': 'a', 'number': 1},
     {'character': 'a', 'number': 2},
     {'character': 'b', 'number': 1},
     {'character': 'b', 'number': 2}]
    """
    return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))

def run_networkclustering (adata, n_neighbors, algorithm) :  
    """
    Runs given algorithm on anndata. Creates new obs column 'algorithm_neighbor'.
    knn = True (default)

    Parameters
    ----------
        n_neighbors (int): range [2,100], default 15
        adata (annData): adata.obs gets extended
        algorithm: see which algorithms work with sc.tl., for example louvain and leiden

    Returns
    -------
        None:  modifies adata
    """
    # if error: file not found then reinstall numpy
    sc.pp.neighbors(adata, n_neighbors = n_neighbors, knn = True) 
    eval("sc.tl." + algorithm + "(adata)")
    adata.obs.rename(columns={algorithm: "{}_{}".format(algorithm, n_neighbors)}, inplace=True)
    

    
def docstringtest() -> int :
    r"""Summarize the function in one line.

    Several sentences providing an extended description. Refer to
    variables using back-ticks, e.g. `var`.

    Parameters
    ----------
    var1 : array_like
        Array_like means all those objects -- lists, nested lists, etc. --
        that can be converted to an array.  We can also refer to
        variables like `var1`.
    var2 : int
        The type above can either refer to an actual Python type
        (e.g. ``int``), or describe the type of the variable in more
        detail, e.g. ``(N,) ndarray`` or ``array_like``.
    *args : iterable
        Other arguments.
    long_var_name : {'hi', 'ho'}, optional
        Choices in brackets, default first when optional.

    Returns
    -------
    type
        Explanation of anonymous return value of type ``type``.
    describe : type
        Explanation of return value named `describe`.
    out : type
        Explanation of `out`.
    type_without_description

    Other Parameters
    ----------------
    only_seldom_used_keyword : int, optional
        Infrequently used parameters can be described under this optional
        section to prevent cluttering the Parameters section.
    **kwargs : dict
        Other infrequently used keyword arguments. Note that all keyword
        arguments appearing after the first parameter specified under the
        Other Parameters section, should also be described under this
        section.

    Raises
    ------
    BadException
        Because you shouldn't have done that.

    See Also
    --------
    numpy.array : Relationship (optional).
    numpy.ndarray : Relationship (optional), which could be fairly long, in
                    which case the line wraps here.
    numpy.dot, numpy.linalg.norm, numpy.eye

    Notes
    -----
    Notes about the implementation algorithm (if needed).

    This can have multiple paragraphs.

    You may include some math:

    .. math:: X(e^{j\omega } ) = x(n)e^{ - j\omega n}

    And even use a Greek symbol like :math:`\omega` inline.

    References
    ----------
    Cite the relevant literature, e.g. [1]_.  You may also cite these
    references in the notes section above.

    .. [1] O. McNoleg, "The integration of GIS, remote sensing,
       expert systems and adaptive co-kriging for environmental habitat
       modelling of the Highland Haggis using object-oriented, fuzzy-logic
       and neural-network techniques," Computers & Geosciences, vol. 22,
       pp. 585-588, 1996.

    Examples
    --------
    These are written in doctest format, and should illustrate how to
    use the function.

    >>> a = [1, 2, 3]
    >>> print([x + 3 for x in a])
    [4, 5, 6]
    >>> print("a\nb")
    a
    b
    """


# code taken from: https://guido.vonrudorff.de/2015/ipython-notebook-code-output-as-markdown/
# code can be found under codesnippets.ipyn as well.


@magics_class
class MarkdownMagics(Magics):

    @cell_magic
    def asmarkdown(self, line, cell):
        buffer = StringIO()
        stdout = sys.stdout
        sys.stdout = buffer
        try:
            exec(cell, globals(), self.shell.user_ns) # changed locals to globals
        except:
            sys.stdout = stdout
            raise
        sys.stdout = stdout
        return HTML("{}".format(markdown(buffer.getvalue(), extensions=['markdown.extensions.extra'])))
           # return buffer.getvalue() + 'test'

#get_ipython().register_magics(MarkdownMagics)