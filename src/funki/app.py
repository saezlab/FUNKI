from dash import Dash
from dash import html
from dash import dcc

from funki import _colors
from pages.home import tab_home
from pages.data import tab_data
from pages.norm import tab_norm
from utils.style import global_style


app = Dash(
    name=__name__,
    title='FUNKI',
)

storage_type = 'session'

app.layout = html.Div(
    children=[
        dcc.Store(id='raw-data', storage_type=storage_type),
        dcc.Store(id='ann-data', storage_type=storage_type),
        dcc.Store(id='proc-data', storage_type=storage_type),
        html.Div(
            html.Img(
                src='assets/logos/funki_logo.svg',
                style={'width': 500, 'padding': 10}
            ),
            style={'width': '100%', 'background-color': _colors['white']}
        ),
        html.Br(),
        dcc.Tabs(
            id='tabs',
            value='tab-home',
            vertical=True,
            children=[
                tab_home,
                tab_data,
                tab_norm,
            ],
            style={
                'padding': 15,
                'width': 150
            },
            colors={
                'background': _colors['white'],
                'border': _colors['teal'],
                'primary': _colors['aqua'],
            }
        )
    ],
    style=global_style,
)

if __name__ == '__main__':
    # Starting the app
    app.run_server(debug=True)