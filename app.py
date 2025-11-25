from dash import Dash, dcc, html
from flame_temp_calculator import Compound_list

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
                options = [{"label": c.name, "value": c} for c in Compound_list],
                value = Compound_list[0]
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

def control_panel(app) -> html.Div:

    return html.Div(
        className = "control-panel",
        children = [
            html.H1(app.title),
            html.Hr(),
            mode_dropdown(),
            html.Div(id = "mode-controls")
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
    app.run()

if __name__ == "__main__":
    main()