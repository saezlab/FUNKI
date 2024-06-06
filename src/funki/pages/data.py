import dash
from dash import html
from dash import dcc

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
            children=[
                'Drag and drop or ',
                html.A('select a file')
            ],
        ),
    ],
    style=tab_style,
    selected_style=tab_selected_style,
)