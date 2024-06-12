import pandas as pd
from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate
from dash.dash_table import DataTable
import plotly.express as px
import plotly.graph_objects as go 
from plotly.subplots import make_subplots

from utils import parse_contents
from utils import serial_to_dataframe
from utils.style import tab_style
from utils.style import tab_selected_style
from utils.style import page_style
from utils.style import header_style


tab_data = dcc.Tab(
    label='Data',
    value='tab-data',
    children=html.Div(
        children=[
            html.H1('Data loading', style=header_style),
            html.Br(),
            html.Div(
                children=[
                    html.Div('Please upload your data file here:'),
                    dcc.Upload(
                        id='upload-data',
                        children=html.Div([
                            'Drag and drop or ',
                            html.A('select a file')
                        ]),
                        style={
                            'width': '80%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '10px'
                        },
                        multiple=False,
                    ),
                    html.Br(),
                    html.Div(
                        children=[
                            dcc.Loading(
                                DataTable(
                                    id='table-data',
                                    fixed_rows={'headers': True, 'data': 0},
                                    fixed_columns={'headers': True, 'data': 1},
                                    style_table={
                                        'maxHeight': 500,
                                        'minWidth': '100%',
                                        'overflowY': 'auto',
                                        'overflowX': 'auto'
                                    },
                                    style_cell={'width': '50px'}
                                )
                            ),
                            dcc.Loading(
                                dcc.Graph(id='plot-data-summary')
                            )
                        ],
                        style={'width': '95%'}
                    ),
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                    'pad': 10
                }
            ),
            html.Div(
                children=[
                    html.Div('Please upload your annotation file here:'),
                    dcc.Upload(
                        id='upload-anndata',
                        children=html.Div([
                            'Drag and drop or ',
                            html.A('select a file')
                        ]),
                        style={
                            'width': '80%',
                            'height': '60px',
                            'lineHeight': '60px',
                            'borderWidth': '1px',
                            'borderStyle': 'dashed',
                            'borderRadius': '5px',
                            'textAlign': 'center',
                            'margin': '10px'
                        },
                        multiple=False,
                    ),
                    html.Br(),
                    html.Div(
                        children=[
                            dcc.Loading(
                                DataTable(
                                    id='table-anndata',
                                    fixed_rows={'headers': True, 'data': 0},
                                    fixed_columns={'headers': True, 'data': 1},
                                    style_table={
                                        'maxHeight': 500,
                                        'minWidth': '100%',
                                        'overflowY': 'auto',
                                        'overflowX': 'auto'
                                    },
                                    style_cell={'width': '50px'}
                                )
                            ),
                            dcc.Loading(
                                dcc.Graph(id='plot-anndata-summary')
                            )
                        ],
                        style={'width': '95%'}
                    ),
                ],
                style={
                    'width': '49%',
                    'display': 'inline-block',
                    'vertical-align': 'top',
                    'pad': 10
                }
            ),
        ],
        style=page_style,
    ),
    style=tab_style,
    selected_style=tab_selected_style,
)

@callback(
    Output('raw-data', 'data'),
    Output('proc-data', 'data'),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    prevent_initial_call=True
)
def load_data(content, filename):
    if filename is None:
        raise PreventUpdate

    df = parse_contents(content, filename)
    serial = {'index': df.index, 'records': df.to_dict('records')}

    return serial, serial

@callback(
    Output('ann-data', 'data'),
    Input('upload-anndata', 'contents'),
    State('upload-anndata', 'filename'),
    prevent_initial_call=True
)
def load_anndata(content, filename):
    if filename is None:
        raise PreventUpdate
    
    df = parse_contents(content, filename)
    serial = {'index': df.index, 'records': df.to_dict('records')}
    
    return serial

@callback(
    Output('table-data', 'columns'),
    Output('table-data', 'data'),
    Output('plot-data-summary', 'figure'),
    Input('raw-data', 'data')
)
def update_table(data):
    if data is None:
        raise PreventUpdate

    df = serial_to_dataframe(data)
    
    fig = px.histogram(df.values.flat, title='Value distribution')
    fig.update_layout(showlegend=False)

    df.reset_index(inplace=True)
    #df = df.head()

    table_columns = [{'name': i, 'id': i} for i in df.columns]
    table_data = df.to_dict('records')

    return table_columns, table_data, fig

@callback(
    Output('table-anndata', 'columns'),
    Output('table-anndata', 'data'),
    Output('plot-anndata-summary', 'figure'),
    Input('ann-data', 'data')
)
def update_anntable(data):
    if data is None:
        raise PreventUpdate
    
    df = serial_to_dataframe(data)

    fig = make_subplots(
        rows=df.shape[1],
        cols=1,
        specs=[[{'type': 'pie'}], [{'type': 'pie'}]]
    )
    
    for i, c in enumerate(df.columns):
        aux = df[c].value_counts()
        pie = go.Pie(labels=aux.index, values=aux.values)
        fig.add_trace(pie, row=i + 1, col=1)

    fig.update_layout(height=250 * df.shape[1])

    df.reset_index(inplace=True)
    #df = df.head()

    table_columns = [{'name': i, 'id': i} for i in df.columns]
    table_data = df.to_dict('records')

    return table_columns, table_data, fig