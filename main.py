# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Main File
# ###################

from dash import dash, dcc, html, Input, Output
from ThermochemicalDataEditor import Compound_list
import plotly.graph_objs as go

##### Setup #####

app = dash.Dash(__name__)

app.layout = html.Div([
    html.Div([
        html.H3("Compound Controls"),

        html.Label("Select Compound"),
        dcc.Dropdown(
            id="compound-dropdown",
            options=[{"label": c.name, "value": c.name} for c in Compound_list],
            value="Methane"
        ),

        html.Label("Select Variable"),
        dcc.Dropdown(
            id="variable-dropdown",
            options=[
                {"label": "Heat of Formation (Hf)", "value": "hf"},
                {"label": "Sensible Heat (Sh)", "value": "sh"},
                {"label": "log Kf", "value": "logKf"}
            ],
            value="hf"
        ),

        html.Label("Display Options"),
        dcc.Checklist(
            id="options-checklist",
            options=[
                {"label": "Show grid", "value": "grid"},
                {"label": "Log scale (y-axis)", "value": "logy"}
            ],
            value=[]
        )
    ],
    style={"width": "25%", "display": "inline-block", "verticalAlign": "top", "padding": "20px"}),

    html.Div([
        dcc.Graph(id="compound-graph", style={"height": "90vh"})
    ],
    style={"width": "70%", "display": "inline-block", "padding": "20px"})
])

##### Callback Logic #####

@app.callback(
    Output("compound-graph", "figure"),
    Input("compound-dropdown", "value"),
    Input("variable-dropdown", "value"),
    Input("options-checklist", "value")
)
def update_graph(selected_compound, selected_variable, options):
    comp = next(c for c in Compound_list if c.name == selected_compound)
    if selected_variable == "hf":
        data = comp.hf_table()
        y_label = "Heat of Formation (kJ/mol)"
    elif selected_variable == "sh":
        data = comp.sh_table()
        y_label = "Sensible Heat"
    else:
        data = comp.logKf_table()
        y_label = "log Kf"

    T, y = data
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=T, y=y, mode="lines+markers", name=selected_compound))

    fig.update_layout(
        title=f"{selected_compound} - {y_label}",
        xaxis_title="Temperature (K)",
        yaxis_title=y_label,
        template="plotly_white",
        showlegend=False
    )

    if "grid" in options:
        fig.update_xaxes(showgrid=True)
        fig.update_yaxes(showgrid=True)
    else:
        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)

    if "logy" in options:
        fig.update_yaxes(type="log")

    return fig

##### Run ######
if __name__ == "__main__":
    app.run(debug=True)