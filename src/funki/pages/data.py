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
from utils import dataframe_to_serial
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
                        id='upload-obs',
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
                                    id='table-obs',
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
                                dcc.Graph(id='plot-obs-summary')
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
    Output('data', 'data', allow_duplicate=True),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    prevent_initial_call=True
)
def load_data(content, filename):
    if filename is None:
        raise PreventUpdate

    serial = dataframe_to_serial(parse_contents(content, filename))

    return serial

@callback(
    Output('data', 'data', allow_duplicate=True),
    Input('upload-obs', 'contents'),
    State('upload-obs', 'filename'),
    State('data', 'data'),
    prevent_initial_call=True
)
def load_obs(content, filename, data):
    if None in (filename, data):
        raise PreventUpdate
    
    serial = dataframe_to_serial(parse_contents(content, filename))
    data.update({'obs': serial})
    
    return data

@callback(
    Output('table-data', 'columns'),
    Output('table-data', 'data'),
    Output('plot-data-summary', 'figure'),
    Input('data', 'data')
)
def update_data_preview(data):
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
    Output('table-obs', 'columns'),
    Output('table-obs', 'data'),
    Output('plot-obs-summary', 'figure'),
    Output('plot-obs-summary' ,'style'),
    Input('data', 'data')
)
def update_obs_preview(data):
    if data is None:
        raise PreventUpdate

    elif 'obs' not in data.keys():
        raise PreventUpdate

    df = serial_to_dataframe(data['obs'])

    fig = make_subplots( # TODO: Improve plotting metadata esp. legend
        rows=df.shape[1],
        cols=1,
        specs=[[{'type': 'pie'}]] * df.shape[1]
    )
    
    for i, c in enumerate(df.columns):
        aux = df[c].value_counts()
        pie = go.Pie(labels=aux.index, values=aux.values)
        fig.add_trace(pie, row=i + 1, col=1)

    height = 250 * df.shape[1]
    fig.update_layout(height=height)

    df.reset_index(inplace=True)
    #df = df.head()

    table_columns = [{'name': i, 'id': i} for i in df.columns]
    table_data = df.to_dict('records')

    return table_columns, table_data, fig, {'height': height}