import streamlit as st
import os, base64, re, sys
from os import makedirs, path
from copy import deepcopy
from IPython.display import display, Markdown   #, Latex # to display Markdown in code chunk
import json, yaml, pickle, uuid                               # numba, logging, random, dill, logging.config
from standard_workflows import *  
sys.path.append('../')
from streamlit_extras.switch_page_button import switch_page
from streamlit_extras.dataframe_explorer import dataframe_explorer 
from st_pages import Page, show_pages, hide_pages
import pandas as pd


#url('https://i.ibb.co/CP9qPhS/FUNKI.png');
##https://i.postimg.cc/vBp8mDdG/funki-human-And-Mice.jpg');#https://i.ibb.co/Hg8yqF0/funki-human-And-Mice.jpg');
 # ":linked_paperclips:"),



class UiVal:
    """Possible values in the UI widgets. (i.e. labels or dropdown values) 
    """
    # 
    PROJ = 'project_id'
    # w_organism
    HUMAN = 'human'
    MOUSE = 'mouse'
    # w_omicstype
    BULKRNA = 'BulkRNAseq'
    PHOSPHO = 'Phosphoproteomics'
    SCRNA = 'SingleCellRNAseq'
    # w_inputformat
    GENES = 'Gene list'
    KINASES = 'Genesymbol_phosphosite list'
    MATRIX = 'Matrix'
    CSV = '.csv'
    H5AD = 'h5ad'
    EXCEL = '.xlsx'
    TSV = '.tsv'

def get_analysis_params(w_organism = UiVal.HUMAN.lower(), w_omicstype = UiVal.BULKRNA.lower()): 
    d = {
    'proj_params': { 
        'proj_id': 'FUNKI',    # always one folder above, projID should be 2 to 8 capital letters or underscores
        'version': 'v01', 
        'paths': {
            'analysis_path': './' ,#path.abspath('../../'),   # path to 'projects' folder or the folder where the proj results shall be saved
            'data_root_path': './data/example_inputs/'#path.abspath('./example_inputs') #+ '/<default>'  # for example path to SDS mounted location: .../mounted/projects/
        },
        'use_pickle_data': False,    # h5ad files are read and then saved as pickle, the pickle files are used from there on
        'datasetname_default': '01',
        'priorKnowledge':{
            'pathways': {
                'progeny':{
                    'top': [300]
                }
            },
            'transcription_factors':{
                'collectri':{
                    'split_complexes': [False]
                }
            },
            'kinase_substrate':{
                'ksn':{}
            }
        },
        'decoupler':{
            'methods': [('ulm', ),],
            #minsize
            #numberOfPermutations
            'meanacts': {
                'groupby': ['vars'],
                'minstd': [0.0]
            }
        },
        'liana': {
            "methods": [["natmi", "connectome", "logfc", "sca", "cellphonedb"]],
            "base": "exp(1)",
            "lig_rec": [["all"]]
        }
    },
    'dataset_params': {
        w_organism: { # organism (human, mouse)
            w_omicstype: { # seqType (scRNA, bulkRNA)
                'priorKnowledge':{
                    'pathways': {
                        'progeny':{
                            'top': [300]
                        }
                    },
                    'transcription_factors':{
                        'collectri':{
                            'split_complexes': [False]
                        }
                    },
                    'kinase_substrate':{
                        'ksn':{
                        }
                    }
                }
            }
        }  
    }
    }
    return d

def set_priorKnwldg(w_omicstype, analysis_params) -> dict:
    """Get priorKnwldg depending on omics type.

    Args:
        w_omicstype (UiVal): widget input
    Return: dict
    """
    def set_priorKnwldg_bykeys(keys:list)-> dict:
        """
        1. Get all priorKnlwdg from project params.
        2. Subset priorKnwldg by given keys. 
        3. Set priorKnwldg for dataset.

        Args:
            keys (List): The priorKnwldg that shall be used.
        Return: dict
        """
        all_priorKnwldg = dict(analysis_params['proj_params']['priorKnowledge'].items())
        subset = sc_funcs.subset_dict(all_priorKnwldg, keys)
        priorKnwldg_dspath = sc_funcs.getpath(analysis_params['dataset_params'], 'priorKnowledge')
        subset = sc_funcs.dict_replace(analysis_params['dataset_params'], subset, priorKnwldg_dspath + ('priorKnowledge',))
        analysis_params['dataset_params'].update(subset)
        return analysis_params

    match w_omicstype:
        case UiVal.BULKRNA | UiVal.SCRNA: 
            return set_priorKnwldg_bykeys(['pathways', 'transcription_factors'])
        case UiVal.PHOSPHO:
            return set_priorKnwldg_bykeys(['kinase_substrate'])

def update_param(param_id, value, dict):
    sc_funcs.dict_replace(dict, value, sc_funcs.getpath(dict, param_id) + (param_id,))


def show_table(data, download_key) -> None:
    """Adds a table and a download button.

    Args:
        data (DataFrame): table
        download_key (String): key
    """
    st.dataframe(data.round(4)) 
    #filtered_df = dataframe_explorer(data.round(4), case=False)
    #st.dataframe(filtered_df, use_container_width=True)
    download_button_str = download_button(data, 'data.csv', 'Download table')
    st.markdown(download_button_str, unsafe_allow_html=True)
    
import re
def sentence_case(sentence):
    """ doesn't work with streamlit... result[:1] is nothing instead of first letter"""
    if sentence != '':
        result = re.sub('([A-Z])', r' \1', sentence)
        result = result[:1].upper() + result[1:].lower()
        return result

button_css = f""" 
        <style>
            #button0 {{
                background-color: rgb(255, 255, 255);
                color: rgb(38, 39, 48);
                padding: 0.25em 0.38em;
                position: relative;
                text-decoration: none;
                border-radius: 4px;
                border-width: 1px;
                border-style: solid;
                border-color: rgb(230, 234, 241);
                border-image: initial;
            }} 
            #button0:hover {{
                border-color: rgb(246, 51, 102);
                color: rgb(246, 51, 102);
            }}
            #button0:active {{
                box-shadow: none;
                background-color: rgb(246, 51, 102);
                color: white;
                }}
        </style> """

#@st.cache_data(experimental_allow_widgets=True) 
def add_results(data, fig, title, download_key) -> None:
    """Add a table with download button and the figure.

    Args:
        data (DataFrame): Data to show
        fig (Plotly): Plot to show
        download_key (String): key
    """
    #title = sentence_case(title)
    st.markdown(f'### {title}')
    col1, col2 = st.columns(2, gap = 'small')
    with col1:
        show_table(data, download_key)     
    with col2:
        if isinstance(fig, str):
            import matplotlib.pyplot as plt
            from PIL import Image
            if fig != 'false':
                with Image.open(fig) as im: 
                    st.image(im)
                    import base64
                    with open(fig, "rb") as image2string:
                        im64 = base64.b64encode(image2string.read())
                with open(fig, 'wb') as im:
                    im64 = im64.decode("utf-8")
                    download_button_str = button_css + f'<a download="image.png" id="button0" href="data:image/png;base64,{im64}">Download</a><br></br>'
                    st.markdown(download_button_str, unsafe_allow_html=True)          
            else:
                st.write('There are too many significant transcription factors. Therefore, no plot is produced.')
        else:
            st.plotly_chart(fig)

    
# make any grid with a function
def make_grid(cols,rows):
    grid = [0]*cols
    for i in range(cols):
        with st.container():
            grid[i] = st.columns(rows)
    return grid



def show_dorothea_res(doro, colNumb = 5, counter=0, mean_act=False):
    if(counter <= len(doro)-1):
        colN = colNumb
        if len(doro)-counter-1 < colN: 
            colN = len(doro)-counter
        cols = st.columns(colN)
        for i in list(range(0,colN)):
                with cols[i]:
                    with open(doro[counter + i], "rb") as pdf_file:
                        base64_pdf = base64.b64encode(pdf_file.read()).decode('utf-8')
                    if(mean_act == True):
                        pdf_display = f'<iframe src="data:application/pdf;base64,{base64_pdf}#view=FitH" width="100%" height="800" type="application/pdf" title="heatmap" data="heatmap.pdf?#zoom=50"></iframe>' 
                    else:
                        pdf_display = f'<embed src="data:application/pdf;base64,{base64_pdf}" width="200" height="170" type="application/pdf">' 
                    title = os.path.basename(doro[counter + i])
                    if(mean_act == True):
                        title = re.search("^(.*)_.*.pdf$", title)[1]
                    else: 
                        title = re.search("^.*_(.*).pdf$", title)[1]
                    st.markdown(f'**{title}**')
                    st.markdown(pdf_display, unsafe_allow_html=True)
        show_dorothea_res(doro, colNumb, counter + colNumb, mean_act)


### FILE UPLOAD ###
def fileupload():
    #import scanpy as sc 
    #data = sc.read("/Users/hanna/Documents/projects/SGUI/CTLA4/v00/analysis/mouse/scRNA/01/data/01.h5ad", cache = True)
    #st.write(data)
    #data

    file = st.file_uploader("Please choose a file", type="h5ad")
    if file is not None:
        #print('huhu')
        #print(file)
        import scanpy as sc
        from io import StringIO
        st.write(file.__dir__())
        st.write(file.name)
        st.write(file.getvalue())
    # Adding a file uploader to accept multiple CSV files
    uploaded_files = st.file_uploader("Please choose a CSV file", accept_multiple_files=True)
    for file in uploaded_files:
        bytes_data = file.read()
        st.write("File uploaded:", file.name)
        st.write(bytes_data)





def print_info(analysis, isSt):
    """ Prints basic information about the analysis. """
    if(isSt):
        cols = st.columns(2)
        with cols[0]:
          st.markdown('**Analysis Parameters**  ', unsafe_allow_html=True)
          st.json(json.dumps(analysis.analysis_params['proj_params']['decoupler'], indent=4, sort_keys=True, default=str))
          st.markdown('**Paths**  ', unsafe_allow_html=True)
          st.json(json.dumps(analysis.get_paths(), indent=4, sort_keys=True, default=str))
        with cols[1]: 
          st.markdown('**All**  ', unsafe_allow_html=True)
          st.json(json.dumps(analysis.analysis_params))
    else:
        display(Markdown('**Analysis Parameters**  '))
        print(json.dumps(analysis.analysis_params['proj_params']['decoupler'], indent=4, sort_keys=True, default=str))
        display(Markdown('**Paths**  '))
        print(json.dumps(analysis.get_paths(), indent=4, sort_keys=True, default=str))

def save_paths(self):
    """ Saves analysis params together with paths to yaml file. """
    self.analysis_params['proj_params']['paths'] = self.get_paths()
    with open(path.join(self.get_paths()['analysis_path'], self.analysis_params["proj_params"]["proj_id"], self.analysis_params["proj_params"]["version"], 'analysis_params_paths.yaml'), 'w+') as file:
        yaml.dump(self.analysis_params, file)
    for data in self.datasets:
        sc_funcs.print_paths(data._paths)
        data.analysis_params['paths'] = data.get_paths()
        with open(path.join(data._paths['analysis_path'], data._paths['datasetpath'], 'analysis_params.yaml'), 'w+') as file:
            yaml.dump(data.analysis_params, file)

def gui_paths(self):
    for data in self.datasets:
        st.text(sc_funcs.print_paths(data._paths))

def read_file(file):
    """ Reads file depending on its extension.

    Args:
        file (UploadedFile): file from st.file_uploader

    Returns:
        data: file content
    """
    ext = os.path.splitext(file.name)[-1].lower()
    match ext:
        case UiVal.CSV:
            return pd.read_csv(file)
        case UiVal.EXCEL:
            return pd.read_excel(file)
        case UiVal.TSV: 
            return pd.read_csv(file, delimiter='\t')

#mygrid = make_grid(5,5)
    #pdfdisp = f'<iframe src="/Abc/drawings/ab-1600.pdf" style="position: absolute; height: 100%; width: 100%;"></iframe>'
    #pdf_display = f'<iframe src="data:application/pdf;base64,{base64_pdf}" width="200" height="200" type="application/pdf"></iframe>'
    #st.markdown(pdf_display, unsafe_allow_html=True)
    #st.markdown(pdfdisp, unsafe_allow_html=True)


def download_button(object_to_download, download_filename, button_text, pickle_it=False):
    """
    ! Code taken from here: https://gist.github.com/chad-m/6be98ed6cf1c4f17d09b7f6e5ca2978f and extended to work for image as well.
    Generates a link to download the given object_to_download.
    Params:
    ------
    object_to_download:  The object to be downloaded.
    download_filename (str): filename and extension of file. e.g. mydata.csv,
    some_txt_output.txt download_link_text (str): Text to display for download
    link.
    button_text (str): Text to display on download button (e.g. 'click here to download file')
    pickle_it (bool): If True, pickle file.
    Returns:
    -------
    (str): the anchor tag to download object_to_download
    Examples:
    --------
    download_link(your_df, 'YOUR_DF.csv', 'Click to download data!')
    download_link(your_str, 'YOUR_STRING.txt', 'Click to download text!')
    """
    import io
    from PIL import Image
    if isinstance(object_to_download, str):#io.BufferedReader):
        #import PIL.Image as Image
        #pil_im = Image.fromarray(object_to_download)
        #b = io.BytesIO()
        #pil_im.save(b, 'png')
        #im_bytes = b.getvalue()
        #st.write(type(im_bytes))
        #img_str = base64.b64encode(im_bytes).decode()
        #st.write(img_str)
        #st.write('hu')
        #st.markdown(get_image_download_link(object_to_download,'image.png','Download'), unsafe_allow_html=True)
        dl_link = custom_css + f'<a download="image.png" id="button0" href="data:image/png;base64,{object_to_download}">Download</a><br></br>'
    else:
        if pickle_it:
            try:
                object_to_download = pickle.dumps(object_to_download)
            except pickle.PicklingError as e:
                st.write(e)
                return None

        else:
            if isinstance(object_to_download, bytes):
                pass

            elif isinstance(object_to_download, pd.DataFrame):
                object_to_download = object_to_download.to_csv(index=False)

            # Try JSON encode for everything else
            else:
                object_to_download = json.dumps(object_to_download)

        try:
            # some strings <-> bytes conversions necessary here
            b64 = base64.b64encode(object_to_download.encode()).decode()

        except AttributeError as e:
            b64 = base64.b64encode(object_to_download).decode()

        button_uuid = str(uuid.uuid4()).replace('-', '')
        button_id = re.sub('\d+', '', button_uuid)

        custom_css = f""" 
            <style>
                #{button_id} {{
                    background-color: rgb(255, 255, 255);
                    color: rgb(38, 39, 48);
                    padding: 0.25em 0.38em;
                    position: relative;
                    text-decoration: none;
                    border-radius: 4px;
                    border-width: 1px;
                    border-style: solid;
                    border-color: rgb(230, 234, 241);
                    border-image: initial;
                }} 
                #{button_id}:hover {{
                    border-color: rgb(246, 51, 102);
                    color: rgb(246, 51, 102);
                }}
                #{button_id}:active {{
                    box-shadow: none;
                    background-color: rgb(246, 51, 102);
                    color: white;
                    }}
            </style> """

        dl_link = custom_css + f'<a download="{download_filename}" id="{button_id}" href="data:file/txt;base64,{b64}">{button_text}</a><br></br>'

    return dl_link


def file_selector(folder_path='.'):
    filenames = os.listdir(folder_path)
    selected_filename = st.selectbox('Select a file', filenames)
    return os.path.join(folder_path, selected_filename)

