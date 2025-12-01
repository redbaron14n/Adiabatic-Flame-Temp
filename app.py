from dash import ALL, Dash, dcc, html, Input, MATCH, Output, State
from domain.compound import Compound
from domain.reaction import Reaction
import plotly.graph_objs as go
from domain.compounds import compounds
from numpy.typing import NDArray

DEFAULT_TEMP: float = 298.15

app = Dash(suppress_callback_exceptions=True)


def create_layout() -> html.Div:

    return html.Div(
        className="app-div",
        style = {"display": "flex", "flexDirection": "row"},
        children=[
            control_panel(),
            html.Div(
                id = "graph-panel-container",
                style = {"width": "75%"},
                children = [
                    graph_panel("compound-graph"),
                    graph_panel("reaction-graph"),
                ]
            ),
        ]
    )


@app.callback(
        Output("compound-graph", "style"),
        Output("reaction-graph", "style"),
        Input("mode-dropdown", "value"),
)
def toggle_graph_visibility(mode: str) -> tuple[dict[str, str], dict[str, str]]:
    if mode == "compound":
        return {"display": "block"}, {"display": "none"}
    else:
        return {"display": "none"}, {"display": "block"}


@app.callback(
    Output("mode-controls-div", "children"),
    Input("mode-dropdown", "value"),
)
def update_mode_controls(mode: str):
    if mode == "compound":
        return compound_controls()
    elif mode == "reaction":
        return reaction_controls()
    else:
        return html.Div("Invalid mode")


def control_panel() -> html.Div:

    return html.Div(
        id="control-panel",
        children=[
            html.H1(app.title),
            html.Hr(),
            mode_dropdown(),
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


def mode_dropdown() -> html.Div:
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


@app.callback(
    Output("compound-graph", "figure"),
    Input("compound-update-graph", "n_clicks"),
    State("compound-selection", "value"),
    State("compound-variable", "value"),
)
def on_compound_graph_update(_, compound_id: str, compound_var: str) -> go.Figure:

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
        case _:
            y_label = ""
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


def graph_panel(graph_id: str) -> html.Div:
    return html.Div(
        className="graph-panel",
        id="graph-area",
        children=[dcc.Graph(id=graph_id, style={"height": "90vh"})],
        style={"width": "70%", "display": "inline-block", "padding": "20px"},
    )


def compound_controls() -> html.Div:

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
            html.Button(id="compound-update-graph", children=["Update Graph"])
        ],
    )


@app.callback(
    Output("reaction-variable", "options"),
    Input("reactant-selection", "value"),
)
def on_reactant_selection(r_ids: list[str]) -> list[dict[str, str]]:
    return [{"label": compounds[r].name, "value": compounds[r].id} for r in r_ids]


@app.callback(
    Output("reactant-ratio-boxes", "children"),
    Input("reactant-selection", "value"),
    Input("reaction-variable", "value"),
)
def update_reactant_ratio_boxes(
    all_reactants: list[str], controlled: str
) -> list[html.Div]:

    boxes = []
    for r in all_reactants:
        if r == controlled:
            continue
        boxes.append(
            html.Div(
                [
                    html.Label(f"{compounds[r].name} Ratio: "),
                    dcc.Input(
                        id={"type": "ratio-input", "compound": r},
                        type="number",
                        min=1,
                        value=1,
                    ),
                ],
                style={"margin-bottom": "8px"},
            )
        )
    return boxes


@app.callback(
    Output("reaction-graph", "figure"),
    Input("reaction-update-graph", "n_clicks"),
    State("reactant-selection", "value"),
    State("reaction-variable", "value"),
    State({"type": "ratio-input", "compound": ALL}, "value"),
)
def on_reaction_graph_update(
    _, r_ids: list[str], controlled: str, ratios: list[float | int]
):
    if not ratios:
        return go.Figure()

    # controlled_reactant: Compound = compounds[controlled]
    # reactants = set(compounds[r] for r in r_ids)
    reactants = set(r for r in r_ids)
    all_reactants = reactants.copy()
    reactants.remove(controlled)
    concentrations = {controlled: 1.0}  # Change to dict[str, float]
    for i, r in enumerate(reactants):
        concentrations[r] = ratios[i]  # Change to dict[str, float]
    temps = {}
    for r in all_reactants:
        temps[r] = DEFAULT_TEMP
    x, t = Reaction(all_reactants, temps).calc_flame_table(
        controlled, concentrations
    )
    y_label = "Flame Temperature (K)"
    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=x, y=t, mode="lines+markers", name="Flame Temperature vs Concentration"
        )
    )
    figure.update_layout(
        title=f"Flame Temperature vs {compounds[controlled].name} Concentration",
        xaxis_title=f"{compounds[controlled].name} Concentration (mol fraction)",
        yaxis_title=y_label,
    )
    return figure


def reaction_controls() -> html.Div:

    return html.Div(
        id="reaction-controls",
        children=[
            html.Label("Select Reactants"),
            dcc.Dropdown(
                id="reactant-selection",
                options=[{"label": c.name, "value": c.id} for c in compounds.values()],
                multi=True,
                value=["Methane", "Oxygen"],
            ),
            html.Label("Select Controlled Reactant"),
            dcc.Dropdown(id="reaction-variable", value="Methane"),
            html.Label("Ratios of Other Reactants"),
            html.Div(id="reactant-ratio-boxes"),
            html.Hr(),
            html.Button(id="reaction-update-graph", children=["Update Graph"])
        ],
    )


if __name__ == "__main__":
    app.title = "Adiabatic Flame Temperature"
    app.layout = create_layout()
    app.run(debug=True)
