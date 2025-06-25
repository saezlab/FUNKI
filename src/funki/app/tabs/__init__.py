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