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


def main() -> None:
    app = Dash()
    app.title = "Adiabatic Flame Temperature"
    app.layout = create_layout(app)
    app.run(debug=True)


def create_layout(app) -> html.Div:

    # return html.Div(className="app-div", children=[control_panel(), graph_panel()])
    return html.Div(className="app-div", children=[control_panel(app)])


def control_panel(app) -> html.Div:

    return html.Div(
        className="control-panel",
        children=[
            html.H1(app.title),
            html.Hr(),
            mode_dropdown(app),
            html.Div(id="mode-controls"),
            html.Hr(),
            html.Button(id="update-graph-button", children=["Update Graph"]),
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


@app.callback(Output("mode-controls", "children"), Input("mode-dropdown", "value"))
def update_mode_controls(mode: str) -> html.Div:

    if mode == "compound":
        return compound_controls()
    elif mode == "reaction":
        return reaction_controls()
    else:
        return html.Div("Invalid mode")


@app.callback(
    Output("compound-selection", "children"), Input("compound-selection", "value")
)
def on_compound_selection(selection: str) -> str:
    _selected_compound = selection
    return selection


def compound_controls() -> html.Div:

    # Define callback for "compound-selection" control here.
    # In callback, set _selected_compound.

    # Define callback for "compound-variable" control here.
    # In callback, set _selected_compound_variable.

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


def reaction_controls() -> html.Div:

    return html.Div(
        id="reaction-controls",
        children=[
            html.Label("Select Reactants"),
            dcc.Dropdown(
                id={"type": "reactant-selection", "index": 0},
                options=[{"label": c.name, "value": c.id} for c in compounds.values()],
                placeholder="Select reactants...",
                multi=True,
            ),
            html.Label("Select Controlled Reactant"),
            dcc.Dropdown(
                id={"type": "controlled-dropdown", "index": 0}, options=[], value=None
            ),
            html.Label("Ratios of Other Reactants"),
            html.Div(id={"type": "reactant-ratios", "index": 0}),
            dcc.Checklist(
                options=[{"label": "Calculate Dissociation", "value": "dissociation"}],
                value=[],
            ),
        ],
    )


@app.callback(
    Output({"type": "controlled-dropdown", "index": MATCH}, "options"),
    Input({"type": "reactant-selection", "index": MATCH}, "value"),
)
def on_reactants_selected(reactants: list[str]) -> list[dict[str, str]]:
    options = []
    if reactants:
        options = [{"label": r, "value": r} for r in reactants]
    return options


@app.callback(
    Output({"type": "reactant-ratios", "index": MATCH}, "children"),
    Input({"type": "reactant-selection", "index": MATCH}, "value"),
    Input({"type": "controlled-dropdown", "index": MATCH}, "value"),
)
def update_ratio_boxes(reactants: list[str], controlled: str) -> list[html.Div]:

    boxes = []
    if not reactants:
        return boxes
    for r in reactants:
        if r == controlled:
            continue
        boxes.append(
            html.Div(
                [
                    html.Label(f"{r} Ratio: "),
                    dcc.Input(
                        id={"type": "reactant-input", "compound": r},
                        type="number",
                        min=1,
                        value=1,
                    ),
                ],
                style={"margin-bottom": "8px"},
            )
        )
    return boxes


def graph_panel() -> html.Div:

    return html.Div(
        className="graph-panel",
        id="graph-area",
        children=[dcc.Graph(id="main-graph", style={"height": "90vh"})],
        style={"width": "70%", "display": "inline-block", "padding": "20px"},
    )


# @app.callback(
#     Output("main-graph", "figure"),
#     Input("mode-dropdown", "value"),
#     # Compound mode states
#     State("compound-selection", "value"),
#     State("compound-variable", "value"),
#     # Reaction mode states
#     State({"type": "reactant-selection", "index": 0}, "value"),
#     State({"type": "controlled-dropdown", "index": 0}, "value"),
#     State({"type": "reactant-input", "compound": ALL}, "value"),
#     State({"type": "reactant-input", "compound": ALL}, "id"),
# )
# def update_graph(
#     mode: str,
#     compound_name: str,
#     variable_name: str,
#     reactants: list[str],
#     controlled: str,
#     ratio_values: list[int],
#     ratio_ids: list[str],
# ) -> go.Figure:
#     if mode == "compound":
#         figure = update_compound_graph(compound_name, variable_name)
#     elif mode == "reaction":
#         figure = update_reaction_graph(reactants, controlled, ratio_values, ratio_ids)
#     else:
#         figure = go.Figure()
#     return figure


def update_compound_graph(compound_name: str, variable_name: str) -> go.Figure:

    compound = compounds[compound_name]
    if variable_name == "hf":
        data = compound.hf_table()
        y_label = "Heat of Formation (kJ/mol)"
    elif variable_name == "sh":
        data = compound.sh_table()
        y_label = "Sensible Heat (kJ/mol)"
    else:
        data = compound.logKf_table()
        y_label = "log Kf"
    t, y = data
    figure = go.Figure()
    figure.add_trace(go.Scatter(x=t, y=y, mode="lines+markers", name=compound_name))
    figure.update_layout(
        title=f"{compound_name} - {y_label}",
        xaxis_title="Temperature (K)",
        yaxis_title=y_label,
        template="plotly_white",
        showlegend=False,
    )
    return figure


# def clean_urg_inputs(
#     reactants: list[str], controlled_var: str, ratios: list[int], ratio_ids: list[str]
# ) -> tuple[set[Compound], Compound, dict[str, int]]:

#     reactant_objs = set()
#     controlled_var_obj = None
#     ratio_dict = {}
#     for c in list(compounds.values()):
#         if c.name in reactants:

#             ###### Something to do with getting dictionaries of ratios of compounds because fuck ass couldnt let me work with that to begin with

#             reactant_objs.add(c)
#         if c.name == controlled_var:
#             controlled_var_obj = c
#     if not controlled_var_obj:
#         controlled_var_obj = list(reactant_objs)[0]
#     return reactant_objs, controlled_var_obj, ratio_dict


# def update_reaction_graph(
#     reactants: list[str],
#     controlled_var: str,
#     ratios: list[int],
#     ratio_ids: list[str],
# ) -> go.Figure:

#     reactant_objs, controlled, ratio_dict = clean_urg_inputs(
#         reactants, controlled_var, ratios, ratio_ids
#     )
#     temp_dict = {}
#     for r in reactant_objs:
#         temp_dict[r] = DEFAULT_TEMP
#     reaction_obj = Reaction(reactant_objs, temp_dict)
#     x, y = reaction_obj.calc_flame_table(controlled, ratio_dict)


if __name__ == "__main__":
    main()
