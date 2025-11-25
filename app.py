from dash import Dash, dcc, html

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

def control_panel(app) -> html.Div:

    return html.Div(
        className = "control-panel",
        children = [
            html.H1(app.title),
            html.Hr(),
            mode_dropdown()
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