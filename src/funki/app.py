import dash
from dash import html
from dash import dcc

from funki import _colors
from pages.home import tab_home
from pages.data import tab_data
from utils.style import global_style


app = dash.Dash(
    name=__name__,
    title='FUNKI'
)

app.layout = html.Div(
    children=[
        dcc.Store(data=None, id='raw-data', storage_type='session'),
        html.Img(
            src='assets/logos/funki_logo.svg',
            style={'width': '50%', 'padding': 10}
        ),
        html.Br(),
        dcc.Tabs(
            id='tabs',
            value='tab-home',
            vertical=True,
            children=[
                tab_home,
                tab_data,
            ],
            style={
                'padding': 15,
            },
            colors={
                'background': _colors['white'],
                'border': _colors['blue'],
                'primary': _colors['aqua'],
            }
        )
    ],
    style=global_style,
)

if __name__ == '__main__':
    # Starting the app
    app.run_server(debug=True)