import dash
from dash import html

from pages.main import main

app = dash.Dash(
    name=__name__,
    title='FUNKI'
)

app.layout = main

if __name__ == '__main__':
    # Starting the app
    app.run_server(debug=True)