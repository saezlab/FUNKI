import numpy as np
import pandas as pd
from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate
from dash.dash_table import DataTable
import decoupler as dc

from utils import serial_to_dataset
from utils import dataset_to_serial
from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style
from funki import _colors
import funki.analysis as fan
import funki.plots as fpl


# ================================== LAYOUT ================================== #

tab_enrichment = dcc.Tab(
    label='Enrichment',
    value='tab-enrichment',
    children=html.Div(
        children=[
            html.H1('Clustering', style=header_style),
            html.Br(),
            html.Div(
                children=[
                    html.Div(
                        children=[
                            'Choose a gene set collection:',
                            html.Br(),
                            dcc.Dropdown(
                                id='gset-collection',
                                options=[
                                    {'label': i, 'value': i}
                                    for i in dc.show_resources()
                                ],
                                searchable=True,
                                clearable=True,
                                style={'width': '80%'},
                            ),
                            html.Br(),
                            dcc.Loading(
                                DataTable(
                                    id='table-gset',
                                    fixed_rows={'headers': True, 'data': 0},
                                    fixed_columns={'headers': True, 'data': 1},
                                    style_table={
                                        'maxHeight': 500,
                                        'minWidth': '100%',
                                        'overflowY': 'auto',
                                        'overflowX': 'auto'
                                    },
                                    style_cell={
                                        'width': 100,
                                        'whiteSpace': 'normal'
                                    }
                                )
                            ),
                            html.Br(),
                            html.Div(
                                id='gset-excl-from-col',
                                hidden=True,
                                children=[
                                    '- Select variable to filter by:',
                                    dcc.Dropdown(
                                        id='gset-excl-from-col-select',
                                        searchable=True,
                                        clearable=True,
                                        style={'width': '80%'},
                                    ),
                                    html.Br(),
                                    html.Div(
                                        id='gset-excl-elems-num',
                                        hidden=True,
                                        children=[
                                            '- Select range of values to keep:',
                                            dcc.RangeSlider(
                                                id='gset-excl-elems-num-select',
                                                min=0,
                                                max=1,
                                                tooltip={
                                                    'always_visible': True,
                                                    'placement': 'top'
                                                },
                                            ),

                                        ]
                                    ),
                                    html.Div(
                                        id='gset-excl-elems-cat',
                                        hidden=True,
                                        children=[
                                            '- Select variable(s) to include:',
                                            dcc.Dropdown(
                                                id='gset-excl-elems-cat-select',
                                                searchable=True,
                                                clearable=True,
                                                multi=True,
                                                style={'width': '80%'},
                                            ),
                                        ]
                                    ),
                                    html.Br(),
                                    html.Button(
                                        'Apply filter',
                                        id='apply-gset-filter',
                                        disabled=True
                                    ),
                                ],
                            ),
                        ],
                        style={
                            'width': '46.5%',
                            'display': 'inline-block',
                            'vertical-align': 'top',
                            'padding-right': 15,
                        }
                    ),
                    html.Div(
                        children=[
                            'Select enrichment method(s):',
                            html.Br(),
                            dcc.Dropdown(
                                id='enrich-methods',
                                options=[
                                    {
                                        'label': r['Name'].rstrip('.'),
                                        'value': r['Function']
                                    }
                                    for i, r in dc.show_methods().iterrows()
                                ],
                                searchable=True,
                                clearable=True,
                                multi=True,
                                style={'width': '90%'},
                            ),
                            html.Br(),
                            '- Select variable to enrich for:',
                            dcc.Dropdown(
                                id='gset-select',
                                searchable=True,
                                clearable=True,
                                style={'width': '80%'},
                            ),
                            html.Br(),
                            html.Button(
                                'Compute enrichment',
                                id='apply-enrichment',
                                disabled=True
                            ),
                            html.Br(),
                            dcc.Loading(dcc.Graph(id='plot-enrich')),
                        ],
                        style={
                            'width': '46.5%',
                            'display': 'inline-block',
                            'vertical-align': 'top',
                            'padding-left': 15,
                        }                    
                    )
                ]
            ),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

# ================================ CALLBACKS ================================= #

@callback(
    Output('table-gset', 'columns'),
    Output('table-gset', 'data', allow_duplicate=True),
    Output('gset-excl-from-col-select', 'options'),
    Output('gset-excl-from-col', 'hidden'),
    Output('gset-select', 'options'),
    Input('gset-collection', 'value'),
    prevent_initial_call=True
)
def load_gset_table(gset):
    if gset is None:
        return None, None, [], True, []

    df = dc.get_resource(gset)

    if len(df) == 0:
        return None, [{0: 'Error downloading the data'}], [], True, []

    table_columns = [{'name': i, 'id': i} for i in df.columns]
    table_data = df.to_dict('records')

    options = [
        {'label': c, 'value': c}
        for c in df.columns
        if c != 'genesymbol'
    ]

    return table_columns, table_data, options, False, options

@callback(
    Output('gset-excl-elems-num-select', 'min'),
    Output('gset-excl-elems-num-select', 'max'),
    Output('gset-excl-elems-num', 'hidden'),
    Output('gset-excl-elems-cat-select', 'options'),
    Output('gset-excl-elems-cat', 'hidden'),
    Output('apply-gset-filter', 'disabled'),
    Input('gset-excl-from-col-select', 'value'),
    Input('table-gset', 'data'),
    prevent_initial_call=True
)
def update_filter(col, gset):    
    df = pd.DataFrame(gset)
    
    # Fallback defaults
    min, max = 0, 1
    options = []
    hid_num, hid_cat = True, True
    dis_button = True

    if not col:
        pass

    elif pd.api.types.is_numeric_dtype(df[col]) and df[col].dtype!= bool:
        max, min = df[col].max(), df[col].min()
        hid_num = False
        hid_cat = True
        dis_button = False

    else:
        options = [
            {'label': str(i), 'value': str(i)}
            for i in df[col].unique()
        ]
        hid_cat = False
        hid_num = True
        dis_button = False

    return min, max, hid_num, options, hid_cat, dis_button

@callback(
    Output('table-gset', 'data', allow_duplicate=True),
    Input('apply-gset-filter', 'n_clicks'),
    State('table-gset', 'data'),
    State('gset-excl-from-col-select', 'value'),
    State('gset-excl-elems-num-select', 'value'),
    State('gset-excl-elems-cat-select', 'value'),
    prevent_initial_call=True
)
def apply_gset_filter(n_clicks, gset, col, rng, cats):
    df = pd.DataFrame(gset)

    if cats:
        userows = df[col].astype(str).isin(cats)
    
    elif rng:
        userows = (df[col] >= rng[0]) & (df[col] <= rng[1])

    else:
        userows = df.index

    return df.loc[userows, :].to_dict('records')

@callback(
    Output('apply-enrichment', 'disabled'),
    Input('enrich-methods', 'value'),
    Input('table-gset', 'data'),
    Input('gset-select', 'value')
)
def update_enrich_button(meth, gset_data, gset):

    return not all([meth, gset_data, gset])

@callback(
    Output('plot-enrich', 'figure'),
    Output('data', 'data', allow_duplicate=True),
    Input('apply-enrichment', 'n_clicks'),
    State('data', 'data'),
    State('table-gset', 'data'),
    State('enrich-methods', 'value'),
    State('gset-select', 'value'),
    prevent_initial_call=True
) # TODO: add human/mouse option
def plot_enrich(n_clicks, data, gset_data, meth, gset):
    if not all([data, gset_data, gset]):
        raise PreventUpdate
    
    net = pd.DataFrame(gset_data)
    dset = serial_to_dataset(data)
    dset.drop_duplicates(subset=['genesymbol', gset], inplace=True)

    fan.enrich(
        dset,
        net,
        methods=[m.lstrip('run_') for m in meth],
        target='genesymbol',
        source=gset,
    )