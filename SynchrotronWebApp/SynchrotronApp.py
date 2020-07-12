import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express       as px
import base64
import pandas as pd
import numpy  as np

desired_indexes = np.arange(0,814,20)
dataframe_path  = "27_446_04_TiAlN_df.pickle"
unpickled_df = pd.read_pickle(dataframe_path).iloc[desired_indexes]
#Convert to degree celcius
unpickled_df['temperature'] = unpickled_df['temperature'] - 273.15


image_filename = 'Icons/materials_chemistry-logo.png' # replace with your own image
encoded_image = base64.b64encode(open(image_filename, 'rb').read())

app = dash.Dash()
app.layout = html.Div(children=[

    html.Div(children=[ 
    html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()), style = {'float': 'right', "width":"120px","padding-right":"40px"}),                        
    html.H1('TiAlN Synchrotron Data', style = {"position":"relative","font-family":"serif","padding-left":"220px","padding-top":"10px"})
    ]),

    html.Div(children=[
    dcc.Slider(id='property-slider', 
               min=0, 
               max=len(desired_indexes)-1, step=1, value=0)
    ], style= {'padding-top': "40px"}),

    html.Div(id='slider-output-container', style = {"padding-left":"20px", "font-size":"30px"}),
    dcc.Graph(id = 'scynchrotron_signal')
])

@app.callback(
    dash.dependencies.Output('slider-output-container', 'children'),
    [dash.dependencies.Input('property-slider', 'value')])
def update_selected_temperature(value):
    desired_index = desired_indexes[int(value)]
    temperature   = unpickled_df.loc[desired_index]['temperature'] 
    return 'In situ Temperature: {:.2f} ℃'.format(temperature)

@app.callback(
    dash.dependencies.Output('scynchrotron_signal', 'figure'),
    [dash.dependencies.Input('property-slider', 'value')])
def update_synchrotron_signal_figure(value):
    desired_index = desired_indexes[int(value)]
    xval   = unpickled_df.loc[desired_index]['q']
    yval   = unpickled_df.loc[desired_index]['intensity'] 

    fig = px.scatter(x=xval, y=yval, template = "simple_white", labels = {"x":"Q", "y":"Intensity (a.u.)"})
    fig.update_xaxes(title_font=dict(size=22))
    fig.update_yaxes(title_font=dict(size=22))
    return fig
if __name__ == '__main__':
    app.run_server(debug=True)