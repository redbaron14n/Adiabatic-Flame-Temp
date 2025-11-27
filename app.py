from dash import Dash, dcc, html, Input, MATCH, Output
from config import Compound_list

def mode_dropdown() -> html.Div:

    return html.Div(
        children = [
            html.Label("Graph Mode"),
            dcc.Dropdown(
                id = "mode-dropdown",
                options = [
                    {"label": "Compound Data", "value": "compound"},
                    {"label": "Reaction Flame Temperature", "value": "reaction"}
                ],
                value = "reaction"
            )
        ]
    )

def compound_controls() -> html.Div:

    return html.Div(
        id = "compound-controls",
        children = [
            html.Label("Select Compound"),
            dcc.Dropdown(
                id = "compound-selection",
                options = [{"label": c.name, "value": c.name} for c in Compound_list],
                value = Compound_list[0].name
            ),
            html.Label("Select Variable"),
            dcc.Dropdown(
                id = "compound-variable",
                options = [
                    {"label": "Constant Pressure Heat Capacity (Cp)", "value": "cp"},
                    {"label": "Standard Entropy (S°)", "value": "s"},
                    {"label": "Total Entropy Change (ΔS)", "value": "ds"},
                    {"label": "Sensible Heat (SH)", "value": "sh"},
                    {"label": "Standard Enthalpy of Formation (ΔHf°)", "value": "hf"},
                    {"label": "Standard Gibbs Free Energy of Formation (ΔGf°)", "value": "gf"},
                    {"label": "log Kf", "value": "logKf"}
                ],
                placeholder = "Select a variable..."
            )
        ]
    )

def reaction_controls(app: Dash) -> html.Div:

    return html.Div(
        id = "reaction-controls",
        children = [
            html.Label("Select Reactants"),
            dcc.Dropdown(
                id = {"type": "reactant-selection", "index": 0},
                options = [{"label": c.name, "value": c.name} for c in Compound_list],
                placeholder = "Select reactants...",
                multi = True
            ),
            html.Label("Select Controlled Reactant"),
            dcc.Dropdown(
                id = {"type": "controlled-dropdown", "index": 0},
                options = [],
                value = None
            ),
            html.Label("Ratios of Other Reactants"),
            html.Div(id = {"type": "reactant-ratios", "index": 0}),
            dcc.Checklist(
                options = [{"label": "Calculate Dissociation", "value": "dissociation"}],
                value = []
            )
        ]
    )

def control_panel(app: Dash) -> html.Div:

    @app.callback(
        Output({"type": "reactant-ratios", "index": MATCH}, "children"),
        Input({"type": "reactant-selection", "index": MATCH}, "value"),
        Input({"type": "controlled-dropdown", "index": MATCH}, "value")
    )
    def update_ratio_boxes(reactants: list[str], controlled: str) -> list[html.Div]:

        boxes = []
        if not reactants:
            return boxes
        for r in reactants:
            if r == controlled:
                continue
            boxes.append(html.Div([
                        html.Label(f"{r} Ratio: "),
                        dcc.Input(
                            id = {"type": "reactant-input", "compound": r},
                            type = "number",
                            min = 1,
                            value = 1
                            )],
                    style = {"margin-bottom": "8px"}))
        return boxes

    @app.callback(
        Output({"type": "controlled-dropdown", "index": MATCH}, "options"),
        Input({"type": "reactant-selection", "index": MATCH}, "value")
    )
    def on_reactants_selected(reactants: list[str]) -> list[dict[str, str]]:
        options = []
        if reactants:
            options = [{"label": r, "value": r} for r in reactants]
        return options

    @app.callback(
        Output("mode-controls", "children"),
        Input("mode-dropdown", "value")
    )
    def update_mode_controls(mode: str) -> html.Div:

        if mode == "compound":
            return compound_controls()
        elif mode == "reaction":
            return reaction_controls(app)
        else:
            return html.Div("Invalid mode")

    return html.Div(
        className = "control-panel",
        children = [
            html.H1(app.title),
            html.Hr(),
            mode_dropdown(),
            html.Div(id = "mode-controls"),
            html.Hr(),
            html.Button(
                id = "update-graph-button",
                children = ["Update Graph"]
            )
        ],
        style = {"width": "25%", "display": "inline-block", "verticalAlign": "top", "padding": "20px"}
    )

def graph_panel() -> html.Div:

    return html.Div(
        className = "graph-panel",
        children = [dcc.Graph(id = "main-graph", style = {"height": "90vh"})],
        style = {"width": "70%", "display": "inline-block", "padding": "20px"}
    )

def create_layout(app: Dash) -> html.Div:

    return html.Div(
        className = "app-div",
        children = [
            control_panel(app),
            graph_panel()
        ]
    )

def main() -> None:

    app = Dash()
    app.title = "Adiabatic Flame Temperature"
    app.layout = create_layout(app)
    app.run(debug=True)

if __name__ == "__main__":
    main()