from dash import ALL, Dash, dcc, html, Input, MATCH, Output, State
from domain.reaction import Reaction
import plotly.graph_objs as go
from domain.compounds import compounds

DEFAULT_TEMP: float = 298.15

# These are the variables for the reaction.
# Several UI inputs set them.
# Another part of the UI shows the result of a reaction based on these inputs.
_selected_compound: str = list(compounds.values())[0].id
_selected_compound_variable: str = ""
_reactants: list = []
# more...


def create_layout(app) -> html.Div:

    return html.Div(
        className="app-div", children=[control_panel(app), graph_panel(app)]
    )


def control_panel(app) -> html.Div:

    @app.callback(
        Output("mode-controls-div", "children"),
        Input("mode-dropdown", "value"),
    )
    def update_mode_controls(mode: str):
        if mode == "compound":
            return compound_controls(app)
        else:
            return html.Div("Invalid mode")

    return html.Div(
        id="control-panel",
        children=[
            html.H1(app.title),
            html.Hr(),
            mode_dropdown(app),
            html.Div(id="mode-controls-div"),
            # html.Hr(),
            # html.Button(id="update-graph-button", children=["Update Graph"]),
        ],
        style={
            "width": "25%",
            "display": "inline-block",
            "verticalAlign": "top",
            "padding": "20px",
        },
    )


def mode_dropdown(app) -> html.Div:
    return html.Div(
        children=[
            html.Label("Graph Mode"),
            dcc.Dropdown(
                id="mode-dropdown",
                options=[
                    {"label": "Compound Data", "value": "compound"},
                    {"label": "Reaction Flame Temperature", "value": "reaction"},
                ],
                value="reaction",
            ),
        ]
    )


def compound_controls(app) -> html.Div:

    @app.callback(Input("compound-selection", "value"))
    def on_compound_selection(selection: str):
        _selected_compound = selection

    @app.callback(Input("compound-variable", "value"))
    def on_compound_variable(cv: str):
        _selected_compound_variable = cv

    return html.Div(
        id="compound-controls",
        children=[
            html.Label("Select Compound"),
            dcc.Dropdown(
                id="compound-selection",
                options=[{"label": c.name, "value": c.id} for c in compounds.values()],
                value=_selected_compound,
            ),
            html.Label("Select Variable"),
            dcc.Dropdown(
                id="compound-variable",
                options=[
                    {"label": "Constant Pressure Heat Capacity (Cp)", "value": "cp"},
                    {"label": "Standard Entropy (S°)", "value": "s"},
                    {"label": "Total Entropy Change (ΔS)", "value": "ds"},
                    {"label": "Sensible Heat (SH)", "value": "sh"},
                    {"label": "Standard Enthalpy of Formation (ΔHf°)", "value": "hf"},
                    {
                        "label": "Standard Gibbs Free Energy of Formation (ΔGf°)",
                        "value": "gf",
                    },
                    {"label": "log Kf", "value": "logKf"},
                ],
                placeholder="Select a variable...",
                value=_selected_compound_variable,
            ),
        ],
    )


def graph_panel(app) -> html.Div:
    return html.Div(
        className="graph-panel",
        id="graph-area",
        children=[dcc.Graph(id="main-graph", style={"height": "90vh"})],
        style={"width": "70%", "display": "inline-block", "padding": "20px"},
    )


if __name__ == "__main__":
    app = Dash()
    app.title = "Adiabatic Flame Temperature"
    app.layout = create_layout(app)
    app.run(debug=True)
