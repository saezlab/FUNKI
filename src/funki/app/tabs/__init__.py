from .home import TabHome
from .data import TabData
from .norm import TabNorm
from .dex import TabDex
from .cluster import TabClust
from .enrich import TabEnrich


TABS = {
    'Home': ('home', TabHome),
    'Data': ('data', TabData),
    'Normalization': ('norm', TabNorm),
    'Clustering': ('clust', TabClust),
    'Diff. Expr.': ('dex', TabDex),
    'Enrichment': ('enrich', TabEnrich),
}

# Stores the corresponding parameter matching between the different attribtues
# of the application tabs and the keys of the DataSet's funki parameters
# dictionary.
# Format is:
# (<tab_name>, <attribute_name>, [sequence of DataSet.uns['funki'] keys])
PARAMS = [
    ('data', 'sample', ['sc_pseudobulk', 'sample_col']),
    ('data', 'group', ['sc_pseudobulk', 'groups_col']),
    ('norm', 'max_genes', ['sc_trans_filter', 'max_genes']),
    ('norm', 'min_genes', ['sc_trans_filter', 'min_genes']),
    ('norm', 'mito_pct', ['sc_trans_filter', 'mito_pct']),
    ('norm', 'size_factor', ['sc_trans_normalize_total', 'target_sum']),
    ('norm', 'log_transform', ['sc_trans_normalize_total', 'log_transform']),
    ('clust', 'harmony', ['harmonize']),
    ('clust', 'harmony_var', ['harmonize', 'vars_use']),
    ('clust', 'embedding_method', ['embedding', 'method']),
    ('clust', 'perplexity', ['embedding', 'perplexity']),
    ('clust', 'min_dist', ['embedding', 'min_dist']),
    ('clust', 'spread', ['embedding', 'spread']),
    ('clust', 'alpha', ['embedding', 'alpha']),
    ('clust', 'gamma', ['embedding', 'gamma']),
    ('clust', 'color_var', ['embedding', 'color_by']),
    ('clust', 'resoultion', ['clustering', 'resolution']),
    ('clust', 'clustering_method', ['clustering', 'alg']),
    ('dex', 'method', ['diff_exp', None, 'method']),
    ('dex', 'contrast_var', ['diff_exp', None, 'design_factor']),
    ('dex', 'groupA', ['diff_exp', None, 'contrast_var']),
    ('dex', 'groupB', ['diff_exp', None, 'ref_var']),
    ('dex', 'thr_logfc', ['diff_exp', None, 'thr_logfc']),
    ('dex', 'thr_pval', ['diff_exp', None, 'thr_pval']),
    ('enrich', 'gset', ['enrich', 'gset']),
    ('enrich', 'org', ['enrich', 'org']),
    ('enrich', 'method', ['enrich', 'method']),
    ('enrich', 'contrast', ['enrich', 'contrast']),
]