import plotly.graph_objects as go
import numpy as np
import pandas as pd 
import dash
from dash import dcc
import dash_bio as dashbio
from dash import html
from dash_bio.utils import xyz_reader
from dash.dependencies import Input, Output, State
from dash import dash_table
from dash.dash_table.Format import Format, Scheme
import dash_bootstrap_components as dbc


def simbapre(kernel, df, index, no):
    #idx = index + 7211
    idx = index + 7094
    kernel_mod = kernel[:,:7093]
    max_idx = np.argpartition(kernel_mod[7094 + index], -no)[-no:]
    return max_idx

def find_idx(compound,df):
    idx = df[df['compound'].str.match(compound)].index.values.astype(int)[0]
    return idx

def data_bars(df, column):
    n_bins = 100
    bounds = [i * (1.0 / n_bins) for i in range(n_bins + 1)]
    col_max = 100
    col_min = -100
    ranges = [
        ((col_max - col_min) * i) + col_min
        for i in bounds
    ]
    styles = []
    for i in range(1, len(bounds)):
        min_bound = ranges[i - 1]
        max_bound = ranges[i]
        min_bound_percentage = bounds[i - 1] * 100
        max_bound_percentage = bounds[i] * 100
        
        style = {
            'if': {
                'filter_query': (
                    '{{{column}}} >= {min_bound}' +
                    (' && {{{column}}} < {max_bound}' if (i < len(bounds) - 1) else '')
                ).format(column=column, min_bound=min_bound, max_bound=max_bound),
                'column_id': column
            },
            'paddingBottom': 2,
            'paddingTop': 2
            }
        if max_bound >= 40:    
            background = (
                """
                    linear-gradient(90deg,
                    white 0%,
                    white 50%,
                    #3D9970 50%,
                    #3D9970 {max_bound_percentage}%,
                    white {max_bound_percentage}%,
                    white 100%)
                """.format(
                    max_bound_percentage=max_bound_percentage,
                    color_above='#3D9970'
                )
            )
        elif ((max_bound < 40) and (max_bound >= 0)) :    
            background = (
                """
                    linear-gradient(90deg,
                    white 0%,
                    white 50%,                    
                    #e9d700 50%,
                    #e9d700 {max_bound_percentage}%,
                    white {max_bound_percentage}%,
                    white 100%)
                """.format(
                    max_bound_percentage=max_bound_percentage,
                    color_above='#e9d700'
                )
            )
        else:
            background = (
                """
                    linear-gradient(90deg,
                    white 0%,
                    white {min_bound_percentage}%,                    
                    #FF4136 {min_bound_percentage}%,
                    #FF4136 50%,
                    white 50%,
                    white 100%)
                """.format(
                    min_bound_percentage=min_bound_percentage,
                    color_above='#FF4136'
                )
            )
        style['background'] = background
        styles.append(style)
        

    return styles



# Declarations
data = pd.read_csv('inh_data.csv')
data['H2evo.220ppm'] = data['H2evo.220ppm'].round(decimals=0)

train = data[data['label'].str.match("train")]
test = data[data['label'].str.match("test")]
new_test = data[data['label'].str.match("new_test")]

untested = data[data['label'].str.match("untested")]
tested = data[data['label'].str.match('|'.join(["train","test","new_test"]))]


## Similarity-based discovery of inhibitors
struc = pd.read_csv('structures_commercial.csv')
struc['IE_krr'] = struc['IE_krr'].round(decimals=0)

sim = np.load('K_c3_g1_commercial.npy')



## Properties
ie = tested['H2evo.220ppm']
mol = xyz_reader.read_xyz(datapath_or_datastring='structures/tris.xyz', is_datafile=True)

sel_idx = find_idx('tris',data)
max_idx = simbapre(sim, data, sel_idx, 5)
struc['sim_val'] = 100 * sim[7094 + sel_idx].round(decimals=3)

sim_struc = struc['Filename'].iloc[max_idx[0]]
mol_sim = xyz_reader.read_xyz(datapath_or_datastring='structures/' + str(sim_struc) + '.xyz', is_datafile=True)

### Scatter plots
fig = go.FigureWidget()
fig.add_trace(go.Scattergl(
    x = untested['CV1'],
    y = untested['CV2'],
    customdata=untested['compound'],
    hovertemplate = '<b>%{customdata}</b>',
    mode='markers',
    marker_symbol = 'cross',
    marker=dict(
        color='lightgray',
        line_width=0,
        size=15
    ),
    name= 'Untested',
    unselected=dict(
        marker=dict(opacity=0.5)
    ),
    selected=dict(
        marker=dict(
            size=25
        )
    )
))

fig.add_trace(go.Scattergl(
    x = tested['CV1'],
    y = tested['CV2'],
    customdata=np.vstack((tested['compound'].values, tested['H2evo.220ppm'].values)).transpose(),
    hovertemplate = '<b>%{customdata[0]}</b><br>%{customdata[1]:.0f} %',
    mode='markers',
    marker=dict(
        color=ie,
        colorscale='BrBG',
        line_width=0,
        size=15,
        colorbar=dict(
            title='IE / %',
            outlinecolor='black',
            outlinewidth=1,
            len=0.4,
            dtick=25,
            x=-0.02,
            yanchor='top',
            y=0.97,
            thickness=20)
    ),
    name= 'Tested',
    unselected=dict(
        marker=dict(opacity=0.5)
    ),
    selected=dict(
        marker=dict(
            size=25
        )
    )
))

fig.update_xaxes(range=[min(tested['CV1'])-0.05, max(tested['CV1'])+0.05])
fig.update_yaxes(range=[min(tested['CV2'])-0.05, max(tested['CV2'])+0.05])


  

### Figure Layout
fig.update_layout(
    autosize=False,
    width=900,
    height=700,
    xaxis = {
        'showgrid':False,
        'zeroline':False,
        'visible': False
        },
    yaxis = {
        'showgrid':False,
        'zeroline':False,
        'visible': False
        },
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    margin=dict(
        l=20,
        r=0,
        b=0,
        t=0,
        pad=4
    ),
    legend=dict(
    yanchor="top",
    y=0.9,
    xanchor="left",
    x=0.1
    ),
    hoverlabel=dict(
        bgcolor="white", 
        font_size=14, 
        font_family="Arial"
    ),
    clickmode='event+select'
)

speck_view_main={
    'resolution': 269,
    'ao': 0.1,
    'outline': 1,
    'atomScale': 0.25,
    'relativeAtomScale': 0.33,
    'bonds': True
}

speck_view_sim={
    'resolution': 269,
    'ao': 0.1,
    'outline': 1,
    'atomScale': 0.25,
    'relativeAtomScale': 0.33,
    'bonds': True
}

speck_legend={
    'resolution': 60,
    #'zoom' : 0.5,
    'ao': 0,
    'outline': 1,
    'atomScale': .9,
    'relativeAtomScale': 1,
    'bonds': True
}

modal = html.Div(
    [
        dbc.Button("Documentation", id="info-open"),
        dbc.Modal(
            [
                dbc.ModalHeader("What is ExChem?"),
                dbc.ModalBody(dcc.Markdown('''
                    ###### About
                    
                    Searching for structure-property relationships is an effective approach to predict yet unknown material properties.
                    *ExChem* allows the exploration of vast areas of chemical space by combining machine learning methods with comprehensive molecular
                    databases.

                    ###### Example: Magnesium Dissolution Modulators<sup>4</sup>
                    
                    Small organic molecules that form complexes with corrosive species accelerating the degradation process have shown great 
                    potential to control the dissolution properties of pure magnesium (Mg) materials and its alloys.<sup>1</sup>
                    However, as the chemical space of small organic molecules is effectively infinite, the most challenging task is to find molecules 
                    with beneficial properties for specific applications.
                    Fortunately, recent studies based on a comprehensive database of magnesium dissolution modulators<sup>1</sup> 
                    revealed that for CPMg220 (commercial purity Mg, containing 220 ppm iron impurities) the molecular structure 
                    correlates well with the corrosion inhibition efficiency (IE).<sup>2,3</sup>

                    To further explore the structure-property relationships in magnesium dissolution modulators, a structure-property landscape was generated
                    using a [SOAP kernel](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5729016/) and [sketch-map](https://sketchmap.org/) for 152 compounds, 
                    of which 78 have already been experimentally tested<sup>1</sup>. On this basis, a kernel ridge regression (KRR) model was trained and used to
                    predict the IEs for a second database of over 7000 commercially available chemicals. Computation of a second SOAP kernel combining both databases
                    allows to screen large areas of chemical space for similar structures based on a selected compound of interest, thus facilitating
                    the search for new corrosion inhibitors.
                    
                    *Note: For some compounds duplicates exist in the reference database. However, due to the training error of the KRR model
                    it may happen that experimental and predicted IEs do not match, although the structures are identical. Generally, we cannot 
                    be held accountable for inaccuratelly predicted values.*

                    ###### How To Use
                    
                    Selection of a point in the sketch-map leads to visualization of the corresponding molecular structure in the dataset. 
                    Additionally, structures similar to the selected are presented in a table 
                    along with their CAS number, similarity value and KRR-predicted IE.  
                    Selecting a table row of interest leads to visualization 
                    of the according molecular structure and its [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) string.
                    Atoms are colored according to the CPK coloring scheme.
                    
                
                    
                    
                    ###### References
                    
                    [[1] Comprehensive screening of Mg corrosion inhibitors, *Corrosion Science* **128** 224–240 (2017)](https://www.sciencedirect.com/science/article/abs/pii/S0010938X17303931)  
                    [[2] Data Science Based Mg Corrosion Engineering, *Frontiers in Materials* **6** 53 (2019)](https://doi.org/10.3389/fmats.2019.00053)   
                    [[3] In silico Screening of Modulators of Magnesium Dissolution, *Corrosion Science* 108245 (2020)](https://doi.org/10.1016/j.corsci.2019.108245)  
                    [[4] Exploring Structure-Property Relationships in Magnesium Dissolution Modulators, *npj Materials Degradation* **5** 2 (2021)](https://www.nature.com/articles/s41529-020-00148-z)
                    ''',dangerously_allow_html=True)
                ),
                dbc.ModalFooter(
                    dbc.Button(
                        "Close", id="info-close", className="ml-auto"
                    )
                ),
            ],
            id="modal-centered",
            centered=True,
            size="xl"
        ),
    ]
)

extra = dbc.Container(fluid=True,children=[dbc.Row(
    [
        dbc.Col(
            dbc.NavItem(
                modal,
                ),style={"list-style-type" : "none","color": "white", "offset" : "20%"}
        ),
        dbc.Col(
            dbc.DropdownMenu(
            children=[
                #dbc.DropdownMenuItem("Examples", header=True),
                dbc.DropdownMenuItem("Magnesium Dissolution Modulators", href="#"),
            ],
            nav=True,
            in_navbar=True,
            label="Datasets",
            toggle_style={"color": "white"},
            style={"list-style-type" : "none"}
            ),
        )        
    ],
    align="center",
    className="g-0"
)],className="g-0",style={'color': 'white', 'width' : '100%' })




navbar = dbc.Navbar(
    dbc.Container(fluid=True,children=
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="https://raw.githubusercontent.com/koerper/ExChem/master/id_logo.png", height="30px"), width="auto"),
                        dbc.Col(dbc.NavbarBrand("ExChem: Explore the Chemical Space", className="ml-2"), width="auto"),
                    ],
                    align="center",
                    className="g-0",
                    style={'width' : '100%'}
                ),
                #href="https://plot.ly",
            ),
            dbc.NavbarToggler(id="navbar-toggler"),
            dbc.Collapse(extra, id="navbar-collapse", navbar=True),
        ]),
        color="dark",
        dark=True,
        #sticky="top",
        fixed="top",        
        style={'color': 'white', 'width' : '100%' }
        )

footer = dbc.Navbar(
    dbc.Container(fluid=True,children=[
        dbc.Row([
            dbc.Col(
                [
                        html.A('Created by '),
                        html.A('Tim Würger'#,href='https://koerper.github.io/',style={"color": "white", "text-decoration": "underline"}
                               ),
                    ], width="auto", lg="auto", xl="auto" 
            ),
            dbc.Col([
                html.A("ExChem: Explore The Chemical Space"),
                #html.A(html.I(className="fab fa-github fa-lg"),href="https://github.com/koerper/ExChem",style={"color": "white", "padding-left" : "3%"}) 
                ], width="auto", lg="auto", xl="auto" 
            ),           
            dbc.Col(
                html.A("© 2020 Tim Würger"), width="auto", lg="auto", xl="auto" 
            )        
        ],  align="center",
            className="g-0",
            justify="between",
            style={'width' : '100%'}
        )
    ]),
        color="dark",
        dark=True,
        #sticky="bottom",
        fixed="bottom",
        style={'color': 'white','height' : '5%', 'width' : '100%' }
)


### App properties
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP,'https://use.fontawesome.com/releases/v5.8.1/css/all.css'])

server = app.server

app.title = 'ExChem'

app.layout = html.Div([
    
    html.Div([
        navbar
    ]),
    html.Div([
        
        html.Div([
            dcc.Markdown('''
                         ##### **Dataset** | Magnesium Dissolution Modulators
                         ''', style={'padding-bottom': '1%'}),
            dcc.Graph(
                id='basic-interactions',
                figure=fig,
                config={
                    'displayModeBar': False
                }
            ),
            html.A(href="https://www.hereon.de",children=[
                    html.Img(src="https://raw.githubusercontent.com/koerper/ExChem/master/assets/hereon_logo.png", height="60px")
                ],
                style={'display': 'inline-block', 'vertical-align' : 'center','padding-top' : '2%','padding-right' : '1%'}
                ),
            html.A(href="https://www.tuhh.de",children=[
                    html.Img(src="https://raw.githubusercontent.com/koerper/ExChem/master/assets/tuhh_logo.png", height="35px")
                ],
                style={'display': 'inline-block', 'vertical-align' : 'center','padding-top' : '2%'}
                ),
            
        ],
        style={'width': '56%', 'display': 'inline-block', 'padding-top' : '1%'}),
        #style={'width': '1000px', 'display': 'inline-block', 'padding-top' : '1%'}),
        html.Div([
            
            html.Div([
                html.H5('Selected structure'),
                html.Div([
                    dashbio.Speck(
                        id='my-speck',
                        data=mol,
                        view=speck_view_main,
                        style={'width' : '269px',
                               'height' : '269px'}
                    )
                ], 
                style={'border': '3px solid','width' : '275px','height' : '275px'}),
                html.P(
                    id='inh_text',
                    children = [],
                    style={'width' : '275px','height' : '60px'}
                )
            ],style={'display': 'inline-block', 'vertical-align' : 'top', 'padding-top' : '3%'}),
                html.Div([
                    html.Img(src="https://raw.githubusercontent.com/koerper/ExChem/master/atoms-legend.png", height="275px")
                ],style={'display': 'inline-block', 'vertical-align' : 'top', 'padding-top' : '8%', 'padding-left' : '54px'}),            
            html.Div([
                html.H5('Proposed similar structure'),                
                html.Div([ 
                    dashbio.Speck(
                        id='sim-speck',
                        data=mol,
                        view=speck_view_sim,
                        style={'width' : '269px',
                               'height' : '269px'}
                    )
                ], 
                style={'border': '3px dashed','width' : '275px','height' : '275px'}),
                html.P(
                    id='sim_text',
                    children = [],
                    style={'width' : '275px','height' : '60px','overflow-wrap': 'break-word'}
                )
            ],style={'display': 'inline-block','vertical-align' : 'top', 'padding-left' : '54px', 'padding-top' : '3%'}),

            html.Div([
                html.Div([
                dbc.Row([
                    dbc.Col([
                       "Top"
                    ],width='auto'),
                    dbc.Col([
                        dcc.Dropdown(
                            id='dropdown',
                            options=[
                                {'label': '5', 'value': '5'},
                                {'label': '10', 'value': '10'},
                                {'label': '20', 'value': '20'},
                                {'label': '50', 'value': '50'}
                            ],
                            value='10',
                            clearable=False
                        )],width={"size": "100px", "offset": "10px"}
                    ),
                    dbc.Col([
                        'similar structures'
                    ],width={"size": "auto", "offset": "10px"}),
                ],align='center',
                className="g-0"
                )
                ],style={'padding-bottom' : '10px'}),
                dash_table.DataTable(
                    id='table',
                    columns=[
                        {"name": 'CAS Number', "id": 'Filename'},
                        {"name": 'Similarity / %', "id": 'sim_val', 'type':'numeric', 'format':Format(precision=3)},
                        {"name": 'Predicted IE / %', "id": 'IE_krr'}
                        ],
                    #selected_rows = max_idx,
                    data=struc.iloc[max_idx].to_dict("records"),
                    selected_rows=[0],
                    style_cell_conditional=[
                        {'if': {'column_id': 'Filename'},
                        'width': '120px',
                        'textAlign': 'left',
                        'minWidth': '120px',
                        'maxWidth': '120px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                         'textAlign': 'center'
                        },
                        {'if': {'column_id': 'IE_krr'},
                        'width': '180px',
                        'minWidth': '180px',
                        'maxWidth': '180px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                        'textAlign': 'center'},
                        {'if': {'column_id': 'sim_val'},
                        'width': '140px',
                        'minWidth': '140px',
                        'maxWidth': '140px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                        'textAlign': 'center'}
                    ],
                    fixed_rows={'headers': True},
                    row_selectable='single',
                    sort_action='native',
                    style_table={
                    'height': '350px', 
                    'width': '700px',
                    'overflowY': 'auto'
                    },
                    css=[{'selector': '.row', 'rule': 'margin: 0'}],
                    style_data_conditional=(
                        data_bars(struc, 'IE_krr')
                    )
                )
            ],style={'vertical-align' : 'top'})
        ],
        #style={'width': '49%', 'display': 'inline-block', 'vertical-align' : 'top', 'padding-left' : '3%'})
        style={'width': '700px', 'display': 'inline-block', 'vertical-align' : 'top'})
    ],style={'font-family' : 'Arial', 'padding-top' : '3%', 'padding-left' : '1%', 'padding-bottom' : '3%'}),
    
    footer
]#,style={'display':'flex', 'flex-direction': 'row'}
)


### Callback functions
@app.callback(
    [Output('my-speck', 'data'), Output('inh_text', 'children')],
    [Input('basic-interactions', 'clickData')]
)
def update_molecule_viewer(clickData):
    if clickData is None:
        mol = xyz_reader.read_xyz(datapath_or_datastring='structures/tris.xyz', is_datafile=True)
        mol_text = ['tris, -72%']
        return mol, mol_text
    point_id = clickData['points'][0]['pointNumber']
    curveNumber = clickData['points'][0]['curveNumber']
    if curveNumber == 0:
        identifier = untested['compound'].iloc[point_id]
        mol_text = [str(identifier)]
    if curveNumber == 1:
        identifier = tested['compound'].iloc[point_id]
        ie_val = tested['H2evo.220ppm'].iloc[point_id]
        mol_text = [str(identifier) + ', ' + str(int(ie_val)) + '%']
    mol = xyz_reader.read_xyz(datapath_or_datastring='structures/' + str(identifier) + '.xyz', is_datafile=True)
    return mol, mol_text




@app.callback(
    [Output('table', 'data'), Output('sim-speck', 'data'), Output('sim_text', 'children')],
    [Input('basic-interactions', 'clickData'),
     Input('table', 'selected_rows'),
     Input('dropdown', 'value')]
)
def update_table(clickData, selected_rows, no_rows):
    if clickData is None:
        sel_idx = find_idx('tris',data)
        max_idx = simbapre(sim, data, sel_idx, int(no_rows))
        
        struc['sim_val'] = 100 * sim[7094 + sel_idx].round(decimals=3)
        table_rows = struc.iloc[max_idx].to_dict("records")
        try:
            sim_struc = struc['Filename'].iloc[max_idx[selected_rows[0]]]
        except:
            selected_rows=[0]
            sim_struc = struc['Filename'].iloc[max_idx[selected_rows[0]]]
        mol_sim = xyz_reader.read_xyz(datapath_or_datastring='structures/' + str(sim_struc) + '.xyz', is_datafile=True)
        sim_text = [struc['Identifier'].iloc[max_idx[selected_rows[0]]]]
        return table_rows, mol_sim, sim_text
    point_id = clickData['points'][0]['pointNumber']
    curveNumber = clickData['points'][0]['curveNumber']
    if curveNumber == 0:
        identifier = untested['compound'].iloc[point_id]
    if curveNumber == 1:
        identifier = tested['compound'].iloc[point_id]
    sel_idx = find_idx(identifier, data)
    
    struc['sim_val'] = 100 * sim[7094 + sel_idx].round(decimals=3)
    max_idx = simbapre(sim, data, sel_idx, int(no_rows))
    table_rows = struc.iloc[max_idx].to_dict("records")
    try:
        sim_struc = struc['Filename'].iloc[max_idx[selected_rows[0]]]
    except:
        selected_rows=[0]
        sim_struc = struc['Filename'].iloc[max_idx[selected_rows[0]]]
    mol_sim = xyz_reader.read_xyz(datapath_or_datastring='structures/' + str(sim_struc) + '.xyz', is_datafile=True)
    sim_text = [struc['Identifier'].iloc[max_idx[selected_rows[0]]]]
    return table_rows, mol_sim, sim_text

@app.callback(
    Output("modal-centered", "is_open"),
    [Input("info-open", "n_clicks"), Input("info-close", "n_clicks")],
    [State("modal-centered", "is_open")],
)
def toggle_modal(n1, n2, is_open):
    if n1 or n2:
        return not is_open
    return is_open

@app.callback(
    Output("navbar-collapse", "is_open"),
    [Input("navbar-toggler", "n_clicks")],
    [State("navbar-collapse", "is_open")],
)
def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

### run server
#if __name__ == '__main__':
#   app.run_server(debug=True, use_reloader=True)  # Turn off reloader if inside Jupyter
