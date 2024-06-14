import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from .analysis import sc_trans_qc_metrics
from .plots import plot_total_counts
from .plots import plot_pct_counts_mito
from .plots import plot_highest_expr
from .plots import plot_counts_vs_n_genes
from .plots import plot_counts_vs_pct_mito
from .plots import plot_n_genes


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