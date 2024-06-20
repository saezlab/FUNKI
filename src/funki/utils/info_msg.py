msg = {
    'upload-data': (
        'Please provide a table where rows are observations (e.g. samples, '
        'cells...) and columns are variables (e.g. genes, proteins...).'
    ),
    'upload-obs': (
        'Please provide a table where rows are observations (e.g. samples, '
        'cells...) and columns are annotation variables (e.g. sample  names, '
        'condition, cell type...).'
    ),
    'filters': (
        'Choosing proper filtering parameters is a highly context-dependent '
        'task. Depending on the experimental set up, these filters can vary a '
        'lot. A proper choice requires some expertise and/or some trial and '
        'error. Therefore, feel free to play with the parameters and observe '
        'the resulting distribution of your data with the help of the plots '
        'below. A good strategy is to start with lax thresholds and move up to '
        'more astringent ones.\n'
        '- Max. genes per cell:\n'
        'Maximum number of genes per cell. Cells with higher number of '
        'measured genes will be excluded from the data set. This filter is '
        'meant to filter out doublets (i.e. when two cells are sequenced as if '
        'they were a single cell).\n'
        '- Min. genes per cell:\n'
        'Minimum number of genes per cell. Cells with a lower number of '
        'measured genes will be excluded from the data set. This filter aims '
        'to filter out ambient RNA detection from empty droplets.\n'
        '- Max. % mito. genes per cell:\n'
        'Maximum percentage of mitochondrial genes allowed per cell. Cells '
        'with a higher percentage of mitochondrial genes will be excluded from '
        'the data set. Expression of mitochondrial genes in the cytoplasm is '
        'generally an indicator of poor sample quality and/or apoptotic or '
        'lysing cells. Nevertheless, the proportion of mitochondrial genes '
        'in a given sample can also depend on the cell type, condition and '
        'experimental procedure. Therefore, these factors should also be taken '
        'into account when choosing a cut-off for this filter.'
    ),
    'norm': (
        'Normalizes each cell by total counts over all genes such that the '
        'total counts of each cell after scaling is equivalent to the '
        'specified size factor. A size factor of one million corresponds to '
        'CPM normalization (counts per million).\n'
        'Although optional, it is usually recommended to also log-transform '
        'the data as it helps reduce skewness and the impact of outliers as '
        'well as bringing the data distribution closer to a normal (Gaussian) '
        'distribution. Specifically, the exact transformation applied to each '
        'value (x) corresponds to: y=log(x+1), where log denotes the natural '
        'logarithm. The addition of one to the original value serves two main '
        'purposes: avoiding division by zero causing the resulting value to be '
        'negative infinite and also, since log(1)=0, non-detected genes remain '
        'valued as zero after the transformation.'
    )
}