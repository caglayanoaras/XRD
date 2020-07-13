import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express       as px
import plotly.graph_objects as go

import base64
import ast
import pandas as pd
import numpy  as np

import XRDLib
import markdown

external_stylesheets = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]
desired_indexes = np.arange(0,814,20)

image_filename = 'Icons/materials_chemistry-logo.png' # replace with your own image
encoded_image = base64.b64encode(open(image_filename, 'rb').read())

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.layout = html.Div(children=[
    html.Div(id='Dataframe', style={'display': 'none'}),

    html.Div(children=[ 
    html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()), style = {'float': 'right', "width":"120px","padding-right":"40px"}),                        
    html.H1('TiAlN Synchrotron Data', style = {"position":"relative","font-family":"serif","padding-left":"220px","padding-top":"10px"})
    ]), 
    
    html.Div(style = {'padding-top': "40px"}, children=[
    dcc.Dropdown(id = "sample_dropdown",
                 options = [{"label":"Sample1", "value":"S1"},{"label":"Sample2", "value":"S2"}], 
                 value ="S1")
    ]),

    html.Div(children=[
    dcc.Slider(id='property-slider', 
               min=0, 
               max=len(desired_indexes)-1, step=1, value=0)
    ], style= {'padding-top': "40px"}),

    html.Div(id='slider-output-container', style = {"padding-left":"20px", "font-size":"30px"}),

    dcc.Graph(id = 'scynchrotron_signal', style = {"height":"60vh"}),

    html.Div([""" Enter peaks:""",html.Div(dcc.Input(id='peak-box', type='text', size = "40")),
              html.Button('Fit', id='button', style = {"height": "40px"}),
              html.Div(id='fit-button', children=[''])
            ], style={"font-size":"20px", "padding-bottom": "10px"}),

    html.Details([
        html.Summary('Fitting Details'),
        html.Div(id = 'detailed_report', children = [])
    ]),

    html.Div(id="graph-container", children = [dcc.Graph(id = 'fit-results')]),

    html.Div(children = [
    html.Div(['''Select plane normal: ''']),
    dcc.RadioItems(id='Select-plane-normal', options=[{'label': i, 'value': i} for i in ['111', '200']], value='111', labelStyle={'display': 'inline-block'})
    ], style = {'text-align':'center', "width":"100%", 'font-size':"20px" } )


])

@app.callback(
    dash.dependencies.Output('Dataframe', 'children'),
    [dash.dependencies.Input('sample_dropdown', 'value')])
def update_sample(value):  
    if value == 'S1':
        dataframe_path = "27_446_04_TiAlN_df.pickle"
    if value == 'S2':
        dataframe_path = "31_447_04_TiAlN_df.pickle"    
    unpickled_df = pd.read_pickle(dataframe_path).iloc[desired_indexes]
    #Convert to degree celcius
    unpickled_df['temperature'] = unpickled_df['temperature'] - 273.15
    return unpickled_df.to_json(date_format='iso', orient='split')


@app.callback(
    dash.dependencies.Output('slider-output-container', 'children'),
    [dash.dependencies.Input('property-slider', 'value'),
     dash.dependencies.Input('Dataframe', 'children')])
def update_selected_temperature(value, df_asjson):
    unpickled_df = pd.read_json(df_asjson, orient = 'split')
    desired_index = desired_indexes[int(value)]
    temperature   = unpickled_df.loc[desired_index]['temperature'] 
    return 'In situ Temperature: {:.2f} ℃'.format(temperature)

@app.callback(
    [dash.dependencies.Output('scynchrotron_signal', 'figure'),#dash.dependencies.Output('fit-button', 'children'),
     dash.dependencies.Output('fit-results', 'figure'),
     dash.dependencies.Output('detailed_report', 'children'),
     dash.dependencies.Output('graph-container', 'style')],
    [dash.dependencies.Input('property-slider', 'value'),
     dash.dependencies.Input('Dataframe', 'children'),
     dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('peak-box', 'value')])
def update_synchrotron_signal_figure(slider_value, df_asjson, n_clicks, box_value):
    unpickled_df = pd.read_json(df_asjson, orient = 'split')
    desired_index = desired_indexes[int(slider_value)]
    xval   = unpickled_df.loc[desired_index]['q']
    yval   = unpickled_df.loc[desired_index]['intensity'] 

    fig = px.scatter(x=xval, y=yval, template = "simple_white", labels = {"x":"Q (1/nm)", "y":"Intensity (a.u.)"})
    fig.update_xaxes(title_font=dict(size=21),tickfont = dict(size=21) )
    fig.update_yaxes(title_font=dict(size=21),tickfont = dict(size=21))
    
    ctx = dash.callback_context

    if not ctx.triggered[0]['prop_id'] == 'button.n_clicks':
        return fig, {}, "" , {'display':'none'}
    else:
        entered_peaks = ast.literal_eval('[' + box_value + ']')
        fit_result = XRDLib.fit_experimental_data(xval, yval,entered_peaks, deg_of_bck_poly = 1, maxfev = 2500)
        fig.add_trace(go.Scatter(x=xval, y=fit_result.best_fit, mode='lines',name='fitted curve',line=dict(color="Crimson"))) 

        fwhm,peak_points,amplitude,fraction, _ = XRDLib.get_profile_data(xval, fit_result) 
        dspacing = [round(2*np.pi/point*10,3) for point in peak_points]
        fit_table = go.Figure(data=[go.Table(header=dict(values=["d_spacing(Å)", 'fwhm', 'peak_points', "amplitude", "fraction"], font_size = 14, height=40),
                              cells=dict(values=[dspacing, fwhm,peak_points,amplitude,fraction], font_size = 14, height=35))
                     ]) 
        return fig, fit_table, [html.P(i) for i in fit_result.fit_report().split('[[Fit Statistics]]')[1].split('\n')], {}

  

if __name__ == '__main__':
    app.run_server(debug=True)