# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Main App File
# ###################

MAX_REACTANTS = 4

import dash
from dash import ALL, dcc, html, Input, Output
from flame_temp_calculator import Compound_list

app = dash.Dash(__name__)

app.layout = html.Div([

    # Mode selection
    html.Div([
        html.Label("Graph Mode"),
        dcc.Dropdown(
            id = "mode-dropdown",
            options = [
                {"label": "Compound Data", "value": "compound"},
                {"label": "Reaction Flame Temperature", "value": "reaction"}
            ],
            value = "reaction"
        )
    ]),

    html.Hr(),

    html.Div(id = "mode-container"),

    html.Hr(),

    dcc.Graph(id = "main-graph")
])

if __name__ == "__main__":
    app.run(debug=True)