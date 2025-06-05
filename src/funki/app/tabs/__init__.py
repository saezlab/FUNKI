from .home import TabHome
from .data import TabData
from .norm import TabNorm
from .dex import TabDex
from .cluster import TabClust
from .enrich import TabEnrich

all_tabs = {
    'Home': ('home', TabHome),
    'Data': ('data', TabData),
    'Normalization': ('norm', TabNorm),
    'Diff. Expr.': ('dex', TabDex),
    'Clustering': ('clust', TabClust),
    'Enrichment': ('enrich', TabEnrich),
}