from dash import ALL, Dash, dcc, html, Input, MATCH, Output, State
from domain.compound import Compound
from domain.reaction import Reaction
import plotly.graph_objs as go
from domain.compounds import compounds
from numpy.typing import NDArray

DEFAULT_TEMP: float = 298.15


def create_layout(app) -> html.Div:

    @app.callback(
        Output("main-graph", "figure"),
        Input("compound-update-graph", "n_clicks"),
        State("compound-selection", "value"),
        State("compound-variable", "value"),
    )
    def on_compound_graph_update(
        test, compound_id: str, compound_var: str
    ) -> go.Figure:

        print("Hello")

        compound: Compound = compounds[compound_id]
        y_vals: NDArray = compound.get_data(compound_var)
        x_vals: NDArray = compound.get_temperatures()
        match compound_var:
            case "cp":
                y_label = "Constant Pressure Heat Capacity (Cp) [J/(mol·K)]"
            case "Hf":
                y_label = "Standard Heat of Formation [kJ/mol]"
            case "SH":
                y_label = "Sensible Heat [kJ/mol]"
            case "logKf":
                y_label = "log Kf"
        figure = go.Figure()
        figure.add_trace(
            go.Scatter(x=x_vals, y=y_vals, mode="lines+markers", name=compound.name)
        )
        figure.update_layout(
            title=f"{compound.name} - {y_label}",
            xaxis_title="Temperature (K)",
            yaxis_title=y_label,
        )
        return figure

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
                value="Methane",
            ),
            html.Label("Select Variable"),
            dcc.Dropdown(
                id="compound-variable",
                options=[
                    # {"label": "Constant Pressure Heat Capacity (Cp)", "value": "cp"},
                    # {"label": "Standard Entropy (S°)", "value": "s"},
                    # {"label": "Total Entropy Change (ΔS)", "value": "ds"},
                    {"label": "Sensible Heat (SH)", "value": "SH"},
                    {"label": "Standard Enthalpy of Formation (ΔHf°)", "value": "Hf"},
                    # {
                    #     "label": "Standard Gibbs Free Energy of Formation (ΔGf°)",
                    #     "value": "gf",
                    # },
                    {"label": "log Kf", "value": "logKf"},
                ],
                placeholder="Select a variable...",
                value="Hf",
            ),
            html.Hr(),
            html.Button(id="compound-update-graph", children=["Update Graph"]),
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
    app = Dash(suppress_callback_exceptions=True)
    app.title = "Adiabatic Flame Temperature"
    app.layout = create_layout(app)
    app.run(debug=True)
