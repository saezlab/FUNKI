import numpy as np
from dash import html
from dash import dcc
from dash import Input
from dash import Output
from dash import State
from dash import callback
from dash.exceptions import PreventUpdate
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
                            'Choose a gene set collection: ',
                            html.Br(),
                            dcc.Dropdown(
                                id='gset-collection',
                                options=[
                                    {'label': i, 'value': i}
                                    for i in dc.show_resources()
                                ],
                                searchable=True,
                                clearable=True,
                            ),
                        ],
                        style={
                            'width': '49%',
                            'display': 'inline-block',
                            'vertical-align': 'top',
                        }
                    ),
                    html.Div(
                        children=[

                        ],
                        style={
                            'width': '49%',
                            'display': 'inline-block',
                            'vertical-align': 'top',
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


