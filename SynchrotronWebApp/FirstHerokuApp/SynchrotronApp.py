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

from WASynchrotron          import Q_to_h3, get_fourier_coefs
from plotly.subplots        import make_subplots
from sklearn.linear_model   import LinearRegression
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
                 options = [{"label":"Sample1 (Claimed to have high defect density)", "value":"S1"},{"label":"Sample2 (Claimed to have low defect density)", "value":"S2"}], 
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
   # html.Div(['''Select plane normal: ''']),
   # dcc.RadioItems(id='Select-plane-normal', options=[{'label': i, 'value': i} for i in ['111', '200']], value='111', labelStyle={'display': 'inline-block'}),
    html.Div([html.Button('Get', id='WAbutton', style = {"height": "40px"}), '''    Warren-Averbach Analysis: '''])
    ], style = {'text-align':'center', "width":"100%", 'font-size':"20px" } ),
    
    html.Div(id = 'WApeakplots-container', children = [dcc.Graph(id = 'WApeakplots111'), dcc.Graph(id = 'WApeakplots200')]),

    html.Div(id = 'WAcoefficientplots-container', children = [html.Div(dcc.Graph(id = 'WAcoefficientplots111'), style={'width':'49%', 'float':'right'}), html.Div(dcc.Graph(id = 'WAcoefficientplots200'), style={'width':'49%', 'float':'left'})], style={'display': 'inline-block'}),

    html.Div(id = 'WAsizeplots-container', children = [html.Div(dcc.Graph(id = 'WAsizeplots111'), style={'width':'49%', 'float':'right'}), html.Div(dcc.Graph(id = 'WAsizeplots200'), style={'width':'49%', 'float':'left'})], style={'height':'40vh', 'display': 'inline-block'}),

    html.Div(id = 'WAstrainplot-container', children = [dcc.Graph(id = 'WAstrainplot')])
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

@app.callback(
    [dash.dependencies.Output('WApeakplots111', 'figure'),
     dash.dependencies.Output('WApeakplots200', 'figure'),
     dash.dependencies.Output('WApeakplots-container', 'style'),
     dash.dependencies.Output('WAcoefficientplots111', 'figure'),
     dash.dependencies.Output('WAcoefficientplots200', 'figure'),
     dash.dependencies.Output('WAcoefficientplots-container', 'style'),
     dash.dependencies.Output('WAsizeplots111', 'figure'),
     dash.dependencies.Output('WAsizeplots200', 'figure'),
     dash.dependencies.Output('WAsizeplots-container', 'style'),
     dash.dependencies.Output('WAstrainplot', 'figure'),
     dash.dependencies.Output('WAstrainplot-container', 'style')],
    [dash.dependencies.Input('WAbutton', 'n_clicks'),
    dash.dependencies.Input('fit-results', 'figure')]
)
def WAanalysis(n_clicks, tabelle):
    ctx = dash.callback_context

    if not ctx.triggered[0]['prop_id'] == 'WAbutton.n_clicks':
        return {},{}, {'display':'none'}, {}, {}, {'display':'none'}, {}, {}, {'display':'none'}, {}, {'display':'none'}
    else:
        if tabelle == {}:
            return {'layout':{'title': "There is no successfull fit!"}},{}, {'display':'none'}, {}, {}, {'display':'none'}, {}, {}, {'display':'none'}, {}, {'display':'none'}
        else:
            def plane_set(plane_family):
                if plane_family == '111':
                    rows = [0,4]
                    legend = ['111', '222']
                if plane_family == '200':
                    rows = [1,5]
                    legend = ['200', '400']
                profile1 = np.array(tabelle['data'][0]['cells']['values'])[:, rows[0]]
                profile2 = np.array(tabelle['data'][0]['cells']['values'])[:, rows[1]]
                d_sp = profile2[0]
                profile2[0] = d_sp*2 
                profile1_spec = Q_to_h3(profile1)
                profile2_spec = Q_to_h3(profile2)

                freq1, coefs1 = get_fourier_coefs(profile1_spec, plot = False)
                freq2, coefs2 = get_fourier_coefs(profile2_spec, plot = False)
                coefslog1 = np.log(coefs1)
                coefslog2 = np.log(coefs2)
                msZn     = []
                msStrain = [0]
                for i,j in zip(coefslog1, coefslog2):
                    msZn.append((i-j)/3/2/(np.pi**2))
                for n,Zn in enumerate(msZn[1:]):
                    msStrain.append(Zn/(n+1)/(n+1))
                msZn[0]   = 0.0 # gets negative value yet so close to zero
                rmsZn     = [np.sqrt(i) for i in msZn]
                rmsStrain = [np.sqrt(i) for i in msStrain]

                Asize = np.exp(np.array(msZn)*2*1*(np.pi**2) +np.log(coefs1))
                reg = LinearRegression().fit( np.arange(5).reshape(-1,1),Asize[:5])   
                dummyx = np.linspace(0,-reg.intercept_/reg.coef_,num=50)
                dummyy = reg.coef_[0]*dummyx + reg.intercept_

                fig_profiles = make_subplots(rows=1, cols = 2)
                fig_profiles.add_trace(go.Scatter(x=profile1_spec.x, y = profile1_spec.y, name=legend[0]), row=1, col=1)
                fig_profiles.add_trace(go.Scatter(x=profile2_spec.x, y = profile2_spec.y, name=legend[1]), row=1, col=2)
                fig_profiles.update_xaxes(title_text="h3", row=1, col=1)
                fig_profiles.update_xaxes(title_text="h3", row=1, col=2)

                fig_coefs = go.Figure()
                count = 0
                for i,j in zip(coefslog1, coefslog2):
                    if count <6:
                        fig_coefs.add_trace(go.Scatter(x=[1,4], y=[i,j], mode= 'lines+markers', line={'dash': 'dash', 'color': 'black'}, name = 'n = {}'.format(count+1)))
                        count +=1
                fig_coefs.update_layout(title_text="Coefficients of {} ".format(plane_family), title_x=0.5)
                fig_coefs.update_xaxes(title_text= "l<sup>2</sup>")
                fig_coefs.update_yaxes(title_text="ln(A<sub>n</sub>)")

                fig_size = go.Figure(layout={'showlegend': False, 'title_text': 'Size coefficients of {}'.format(plane_family), 'title_x':0.5, 'height':500})
                fig_size.add_trace(go.Scatter(x=list(range(200)), y=Asize[:200], mode= 'lines+markers', line={'dash': 'solid', 'color': 'black'}))
                fig_size.add_trace(go.Scatter(x=[0, -reg.intercept_/reg.coef_[0]], y=[1, 0], mode= 'lines+markers', line={'dash': 'dash', 'color': 'red'}))
                fig_size.update_xaxes(title_text='n')
                fig_size.update_yaxes(title_text='A<sub>s</sub>')
                return fig_profiles, fig_coefs, fig_size, rmsStrain

            fig_profiles111, fig_coefs111, fig_size111, rmsStrain111 = plane_set('111')
            fig_profiles200, fig_coefs200, fig_size200, rmsStrain200 = plane_set('200')

            fig_strain = go.Figure(layout  = {'title_text': 'Microstrain'})
            fig_strain.update_xaxes(title_text= "L (Å)")
            fig_strain.update_yaxes(title_text="RMS Strain (%)")
            d111 = np.array(tabelle['data'][0]['cells']['values'])[0][0]
            d200 = np.array(tabelle['data'][0]['cells']['values'])[0][1]
            LL200 = [i*d200 for i in range(17)]
            LL111 = [i*d111 for i in range(15)]
            fig_strain.add_trace(go.Scatter(x=LL200[1:], y= [i*100 for i in rmsStrain200[1:17]], name ='200' ))
            fig_strain.add_trace(go.Scatter(x=LL111[1:], y= [i*100 for i in rmsStrain111[1:15]], name ='111' ))
            
            return fig_profiles111,fig_profiles200,{}, fig_coefs111, fig_coefs200, {}, fig_size111, fig_size200, {}, fig_strain, {}

if __name__ == '__main__':
    app.run_server