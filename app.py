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

@app.callback(
    Output("mode-container", "children"),
    Input("mode-dropdown", "value")
)

def display_mode(mode):

    if mode == "compound":
        return compound_controls()
    elif mode == "reaction":
        return reaction_controls()
    
def compound_controls():

    return html.Div([
        html.Label("Compound"),
        dcc.Dropdown(
            id="compound-dropdown",
            options = [{"label": c.name, "value": c.name} for c in Compound_list],
            value = Compound_list[0].name
        ),

        html.Label("Variable"),
        dcc.Dropdown(
            id = "variable-dropdown",
            options = [
                {"label": "Constant Pressure Specific Heat (Cp)", "value": "cp"},
                {"label": "Standard Entropy (S°)", "value": "s"},
                {"label": "Total Entropy Change (ΔS)", "value": "ds"},
                {"label": "Sensible Heat (SH)", "value": "sh"},
                {"label": "Heat of Formation (ΔHf)", "value": "hf"},
                {"label": "Gibbs Free Energy (ΔGf)", "value": "gf"},
                {"label": "log Kf", "value": "logKf"}
            ],
            value = "sh"
        )
    ])

def reaction_controls():

    return html.Div([
        html.Label("Number of Reactants"),
        dcc.Dropdown(
            id = "num-reactants",
            options = [{"label": str(i), "value": i} for i in range(1, MAX_REACTANTS+1)],
            value = 2
        ),

        html.Div(id = "reactant-selectors"),

        html.Label("Reactant with Variable Concentration"),
        dcc.Dropdown(
            id = "variable-reactant"
        ),

        html.Div(id="fixed-reactant-inputs"),

        dcc.Checklist(
            id = "dissociation-check",
            options = [{"label": "Enable dissociation", "value": "diss"}],
            value = []
        )
    ])

@app.callback(
    Output("reactant-selectors", "children"),
    Input("num-reactants", "value")
)

def build_reactant_selectors(n: int):

    return [
        html.Div([
            html.Label(f"Reactant {i+1}"),
            dcc.Dropdown(
                id = {"type": "reactant-name", "index": i},
                options = [{"label": c.name, "value": c.name} for c in Compound_list],
            )
        ])
        for i in range(n)
    ]

if __name__ == "__main__":
    app.run(debug=True)