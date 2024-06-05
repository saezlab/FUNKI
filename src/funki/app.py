import dash
from dash import html
from dash import dcc

from pages.home import tab_home
from pages.data_load import tab_data_load

app = dash.Dash(
    name=__name__,
    title='FUNKI'
)

app.layout = html.Div(
    children=[
        html.Img(
            src='assets/logos/funki_logo.svg',
            style={'width': '50%', 'padding': '10px'}
        ),
        html.Br(),
        dcc.Tabs(
            id='tabs',
            value='tab-home',
            vertical=True,
            children=[
                tab_home,
                tab_data_load,
            ]
        )
    ]
)

if __name__ == '__main__':
    # Starting the app
    app.run_server(debug=True)