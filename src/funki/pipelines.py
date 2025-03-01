import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .analysis import sc_trans_qc_metrics
from .analysis import diff_exp
from .analysis import enrich
from .plots import plot_total_counts
from .plots import plot_pct_counts_mito
from .plots import plot_highest_expr
from .plots import plot_counts_vs_n_genes
from .plots import plot_counts_vs_pct_mito
from .plots import plot_n_genes
from .plots import plot_dex
from .plots import plot_enrich


def sc_quality_control(data):
    '''
    Computes QC metrics on a single-cell data set and generates several plots
    to visualize them. Generates a multipanel figure with the follwoing plots:

    * Box plot with highest expression genes
    * Violin plot with number of genes per cell
    * Violin plot with total counts per gene
    * Violin plot with the percentage of mitochondrial genes per cell
    * Scatter plot of total counts vs. percentage of mitochondrial genes
    * Scatter plot of total counts vs. number of genes

    :param data: The data set from which to compute the QC metrics
    :type data: :class:`funki.input.DataSet`

    :returns: The figure contataining the resulting plot with multiple panels
        for different metrics and comparisons
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    data.var['mito'] = data.var_names.str.startswith('MT-')
    data = sc_trans_qc_metrics(data)

    fig = make_subplots(rows=2, cols=3, subplot_titles=[' '] * 6)

    plt1 = plot_highest_expr(data, top=5)
    plt2 = plot_n_genes(data)
    plt3 = plot_total_counts(data)
    plt4 = plot_pct_counts_mito(data)
    plt5 = plot_counts_vs_pct_mito(data)
    plt6 = plot_counts_vs_n_genes(data)

    for i, plt in enumerate([plt1, plt2, plt3, plt4, plt5, plt6]):
        c, r = (i % 3) + 1, (i // 3) + 1
        fig.add_trace(plt.data[0], row=r, col=c)
        fig.update_xaxes(
            title_text=plt.layout['xaxis']['title']['text'],
            row=r,
            col=c
        )
        fig.update_yaxes(
            title_text=plt.layout['yaxis']['title']['text'],
            row=r,
            col=c
        )
        fig.layout.annotations[i].update(text=plt.layout['title']['text'])

    return fig

def differential_expression(
    data,
    design_factor,
    contrast_var,
    ref_var,
    logfc_thr=1.0,
    fdr_thr=0.05,
    n_cpus=8
):
    '''
    Computes differential expression analysis on the provided data based on a
    given design factor and both the contrast and reference variables (e.g.
    treatment and control). Then generates the resulting volcano plot based on
    the desired thresholds.

    :param data: The data from which to compute the differential expression
    :type data: :class:`funki.input.DataSet`
    :param design_factor: Name of the column containing the variables which the
        contrasting samples are assigned. The column must be present in the
        ``data.obs`` table
    :type design_factor: str
    :param contrast_var: The variable value(s) that defines the samples that are
        to be contrasted against the reference (e.g. ``'treatment'``). The value
        must be present in the specified ``design_factor`` column
    :type contrast_var: any | list[any]
    :param ref_var: The variable value(s) that defines the refence samples (e.g.
        ``'control'``). The value must be present in the specified
        ``design_factor`` column
    :type ref_var: any | list[any]
    :param logfc_thr: Threshold for signifacnce based on the log2(FC) value,
        defaults to ``1.0``
    :type logfc_thr: float, optional
    :param fdr_thr: Threshold for signifacnce based on the FDR value, defaults
        to ``0.05``
    :type fdr_thr: float, optional
    :param n_cpus: Number of CPUs used for the calculation, defaults to ``8``
    :type n_cpus: int, optional

    :returns: The figure contataining the resulting scatter plot
    :rtype: `plotly.graph_objs.Figure`_

    .. _plotly.graph_objs.Figure: https://plotly.com/python-api-reference/gener\
        ated/plotly.graph_objects.Figure.html
    '''

    diff_exp(data, design_factor, contrast_var, ref_var, n_cpus=n_cpus)
    
    return plot_dex(data, logfc_thr=logfc_thr, fdr_thr=fdr_thr)

def enrichment_analysis(
    data,
    net,
    methods=None,
    source=None,
    target=None,
    weight=None,
    top=10,
    **kwargs
):
    '''
    Performs enrichment analysis using `Decoupler`_ based on a given network
    (e.g. gene set collection) and statistical method(s) and returns a
    figure with the consensus score across methods for the enrichment results.

    :param data: The data set from which to perform the enrichment
    :type data: :class:`funki.input.DataSet`
    :param net: The network linking the features of the data to the attributes
        (e.g. pathways, gene sets, transcription factors, etc.)
    :type net: `pandas.DataFrame`_
    :param methods: Which statistical method(s) to use in order to compute the
        enrichment, defaults to ``None``. If none is provided, uses ``'mlm'``,
        ``'ulm'`` and ``'wsum'``. The option ``'all'`` performs all methods. To
        see all the available methods, you can run `decoupler.show_methods()`_
        function
    :type methods: NoneType | str | list[str]
    :param source: Column name from the provided ``net`` containing the gene
        sets to enrich for.
    :type source: str
    :param target: Column name from the provided ``net`` containing the gene set
        components (e.g. gene/protein names) that can be mapped back to the data
        set variable names.
    :type target: str
    :param weight: Defines the column in the network containing the weights to
        use in the enrichment, defaults to ``None``.
    :type weight: NoneType | str
    :param top: Number of top enriched gene sets to display based on their
        consensus score. If a negative number is provided, the bottom ones will
        be displayed instead.
    :type top: int
    :param \*\*kwargs: Other keyword arguments that passed to
        `decoupler.decouple()`_ function
    :type \*\*kwargs: optional

    :returns: ``None``, results are stored inplace of the passed ``data``
        object, which is a :class:`funki.input.DataSet` instance. Estimates,
        p-values and consensus scores (in case of multiple methods) are stored
        as part of the ``obsm`` attribute of the object.

    :rtype: NoneType

    .. _Decoupler: https://decoupler-py.readthedocs.io/en/latest/index.html
    .. _pandas.DataFrame: https://pandas.pydata.org/docs/reference/api/pandas.D\
        ataFrame.html
    .. _decoupler.show_methods(): https://decoupler-py.readthedocs.io/en/latest\
        /generated/decoupler.show_methods.html#decoupler.show_methods
    .. _decoupler.decouple(): https://decoupler-py.readthedocs.io/en/latest/gen\
        erated/decoupler.decouple.html#decoupler.decouple
    '''
        
    enrich(
        data,
        net,
        methods=methods,
        source=source,
        target=target,
        weight=weight,
        **kwargs
    )
    return plot_enrich(data, top=top)