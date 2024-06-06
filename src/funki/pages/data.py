import dash
from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate

from utils.style import tab_style
from utils.style import tab_selected_style
import funki

tab_data = dcc.Tab(
    label='Data',
    value='tab-data',
    children=[
        html.H1('Data loading'),
        html.Div('Please upload your data file here:'),
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and drop or ',
                html.A('select a file')
            ]),
            style={
                'width': '100%',
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
    ],
    style=tab_style,
    selected_style=tab_selected_style,
)

@callback(
    Output('raw-data', 'data'),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename'),
    prevent_initial_call=True
)
def load_data(contents, filename):
    if filename is None:
        raise PreventUpdate

    return funki.input.read(filename)
