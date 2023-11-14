import subprocess 
#import webbrowser, glob, os, subprocess              #json, yaml, base64, re, sys, numba, logging, random, dill, logging.config
from os import makedirs, path       
#from copy import deepcopy
import streamlit as st
st.set_page_config(layout="wide")
# funki modules
from standard_workflows import * 
#from gui import pages, utils

from gui.utils import utilities as util
from gui.utils import bulkRNA_utils as bulk
from gui.utils import web_utils as web
from gui.utils import analysispage_utils as analysispage
# Logging
import logging
from io import StringIO

##########################################################################################################################################
# st.success(f'The analysis will be run for **{w_organism}** with **{w_omicstype}** data in the format **{w_inputformat}**.')
# #st.session_state.genelist = "Please wait while your data is processing..."
# sparse matrix wrong datatype object instead of float in df coming from pd.read_csv
# row t, cols genes -> row t, cols pathways
# #d.values.astype(float).index

# Code for debugging: 
# data = pd.read_csv('../../example_inputs/differential_stats.csv')
# ap = util.get_analysis_params()
# priorKnwldg = ap['dataset_params']['human']['bulk rnaseq']['priorKnowledge'].keys()
# for resource in priorKnwldg:
#     l = sc_funcs.deep_get(ap, sc_funcs.getpath(ap, resource, search_value = False))
#     print(l)

# counts = csr_matrix(d)
# adata = ad.AnnData(counts)
# adata
# adata.obs_names = [f"Cell_{i:d}" for i in d.index] 
# adata.var_names = [f"Gene_{i:d}" for i in d.columns] 
# print(adata.obs_names[:10])


# r = result_estimate.transpose()
# r = r.head(1)
# adata = ad.AnnData(r, obs = pd.DataFrame(index = r.index), var=pd.DataFrame(index=r.columns), dtype=np.float32)

# adata = ad.AnnData(d, obs=pd.DataFrame(index=d.index), var=pd.DataFrame(index=d.columns))
# TODO: 
# - REQ1: Ora background genes: we only need the number of background genes and not the actual genes but it's more intuitive for the people if they can just upload the genelist instead of counting the genes.
# - Same project and datasetname -> read in analysis params file if exists and execute analsis accordingly on 'submit'
# - Add new dataset
# - Choose already processed datasets that shall show their results in the tabs
# - keep plotting function very general to use it with varous input data (left col for table, right col for plot)
# - if not 'Run' button, don't change any params or results.
# - Extend to multiple datasets at once? Not really necessary. Maybe one tab per ds results?
###################################################################################################################
if 'isScEnabled' not in st.session_state:
    st.session_state.isScEnabled = False

def show_advanced_options(ap, organism, omicstype):
    with st.form('advanced_params', clear_on_submit=False):
        st.caption('Advanced Settings')
        st.warning("This feature is work in progress. You can see changes, that you do here, in the 'Analysis Parameters' tab but ony single cell data gets saved so far.")
        w_datasetname_default = ap['proj_params']['datasetname_default']
        w_projname_default = ap['proj_params']['proj_id']
        w_projname    = st.text_input(UiVal.PROJ, placeholder = w_projname_default)
        w_datasetname = st.text_input('Dataset Name', placeholder = w_datasetname_default)  #datasetnames)
        w_resultspath = st.text_input('Results Path', placeholder = '/Users/MaxMustermann/Documents/')
        w_topn = st.text_input('top n', placeholder = '300, 500')
        w_save_parameters = st.form_submit_button('Save Parameters')
        
        if(w_save_parameters):
            if w_projname != '':
                util.update_param('proj_id', w_projname, ap)
            else:
                w_projname = w_projname_default

            if w_datasetname == '':
                w_datasetname = w_datasetname_default 
            ap['dataset_params'][w_organism][w_omicstype][w_datasetname] = {}

            if w_topn != '':
                ap['dataset_params'][organism][omicstype]['priorKnowledge']['pathways']['progeny'] = {'top': [int(x) for x in list(w_topn.split(','))]}

            ap['proj_params']['paths']['analysis_path'] = w_resultspath
            if w_resultspath != '':
                st.success(f'The results will be saved in **{w_resultspath}/{w_projname}/{w_datasetname}/**.')
                st.session_state.isScEnabled = True
            else:
                st.warning('Please provide a "Results Path" if the results shall be saved.')        
                
        st.session_state.ap = ap # as long as the form was not sent, 'ap' doesn't change


def fill_tab4(ap):
    ##################    
    ### Parameters ###
    ##################
    with tab4:
        import pandas as pd
        st.warning('This feature is work in progress.')

        # Extract params
        st.write('### Chosen Analysis Parameters')

        organism = list(ap['dataset_params'].keys())[0]
        omicstype = list(ap['dataset_params'][organism].keys())[0]
        datasetnames = [k for k in list(st.session_state.ap['dataset_params'][w_organism][w_omicstype].keys()) if k != 'priorKnowledge']
        if len(datasetnames) >= 1:
            datasetnames = datasetnames[0]
        else: 
            datasetnames = ''
        analysispath = ap['proj_params']['paths']['analysis_path']
        priorKnowledge = ap['dataset_params'][organism][omicstype]['priorKnowledge']
        proj_id = ap['proj_params']['proj_id']
        
        data = {'ProjectID': proj_id,
                 'Organism': organism,
                 'OmicsType': omicstype, 
                 'DatasetID': datasetnames,
                 'AnalysisPath': analysispath}
        
        st.dataframe(data)
        data = {
                 'PriorKnowledge': priorKnowledge}
        st.dataframe(data)
        
        st.write('### Internal Representation of Parameters')
        st.write(ap)
        st.write('### Merged')
        if(datasetnames != ''):
            st.write(sc_funcs.merge_dicts(ap['proj_params'], ap['dataset_params'][organism][omicstype][datasetnames]) )
        st.write('### Session State')
        st.write(st.session_state)
        st.write('### Paths')
        #st.markdown(sc_funcs.print_paths(st.session_state.ap['proj_params']['paths']))



        
#--------LOGGING--------#
logger = logging.getLogger(__name__)
#logger.handlers.clear()
print("handlers", logger.handlers)
if not len(logger.handlers):
    log_capture_string = StringIO()
    fhandler = logging.FileHandler(filename='funki.log', mode='a')
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fhandler.setFormatter(formatter)
    fhandler.name = "funkifile_handler"
    logger.addHandler(fhandler)
    logger.setLevel(logging.INFO)


#--------INIT--------#
UiVal = util.UiVal
web.init_page()
analysispage.init_page()

if 'analysiscount' not in st.session_state:
    st.session_state.analysiscount = 0

tab1, tab2, tab3, tab4 = st.tabs(['Analysis', 'Results 1st dataset', 'Results 2nd dataset', 'Parameter Choices'])
#tabs = st.tabs(['Analysis'] + [f'analysis{x}' for x in range(st.session_state.analysiscount)] + ['ParameterChoices'])

################
### Analysis ###
################
with tab1:
    st.caption('Describe your data')
    paramcol1, paramcol2 = st.columns(2)

    #---- DATASET SPECIFIC PARAMS ----#
    with paramcol1:
#w_datasetname_default = '01'
        w_organism    = st.selectbox('Organism :bust_in_silhouette: :mouse2:', (UiVal.HUMAN, UiVal.MOUSE))
        w_omicstype   = st.selectbox('Omics Type', (UiVal.BULKRNA, UiVal.PHOSPHO, UiVal.SCRNA,)) 
        # adjust wording (gene/kinase)
        match w_omicstype: 
            case UiVal.BULKRNA:
                inputformats = (UiVal.GENES, UiVal.MATRIX)
            case UiVal.PHOSPHO:
                inputformats = (UiVal.KINASES, UiVal.MATRIX)
            case UiVal.SCRNA:
                inputformats = (UiVal.H5AD,)

        w_inputformat = st.selectbox('Input format  :page_with_curl:', inputformats)
        w_testdata    = st.checkbox('Use test data (human) :bar_chart:.\n\n*(Check the needed input formats*)')

        if w_omicstype == UiVal.SCRNA:
            st.warning("To use this option you need a .h5ad file where the adata raw field contains the log transformed data. Furthermore, you need to fill in the 'Input Path' and 'Results Path'. The results will be saved directly to the given results path. ")


    #---- INIT ANALYSIS_PARAMS ----#
#if 'ap' not in st.session_state:
    analysis_params = util.get_analysis_params(w_organism, w_omicstype)
    util.set_priorKnwldg(w_omicstype, analysis_params)


    #---- ADVANCED OPTIONS ----#        
    with paramcol2:
        # Project Specific Params
        w_show_all_opts = st.checkbox('Show all options')
        if w_omicstype == UiVal.SCRNA:
            w_show_all_opts = True
        if(w_show_all_opts): 
            show_advanced_options(analysis_params, w_organism, w_omicstype)
        else:
            st.session_state.ap = analysis_params
    #---- Get/Prepare Data ----#

    datasets = list()
    #aps = {'organism': list(st.session_state.ap['dataset_params'].keys())[0]}
    #aps.update({'omicstype': list(st.session_state.ap['dataset_params'][aps['organism']].keys())[0]}) # analysis params    
    if not w_testdata:
        if not (w_organism == UiVal.MOUSE) & (w_omicstype == UiVal.SCRNA):
            if not (w_omicstype == UiVal.SCRNA) & (st.session_state.isScEnabled == False):
                datasets = bulk.get_data(w_inputformat)
                if len(datasets) != 0:       
                    def get_acts_perDs(datasets):
                        with tab2:
                            cols = st.columns(len(datasets)) 
                            for i in range(0, len(datasets)):
                                with cols[i]:
                                    bulk.get_acts(datasets[i])
                    get_acts_perDs(datasets)
        else:
            st.warning("You can't analyse single cell mouse data with FUNKI at the moment. If you want to do so, please write to Hanna and she'll notify you when this problem is fixed")
    else: # let gettestdata return 'datasets'
        datasets = bulk.get_testdata(w_inputformat, w_omicstype, datarootpath = st.session_state.ap['proj_params']['paths']['data_root_path'])
        #st.write('The following data will be used for the analysis: ')
        st.write(datasets['datasetname'])
        #st.write(datasets['data'])
        bulk.get_acts(datasets)
    
fill_tab4(st.session_state.ap)
    




################
### Datasets ###
################
#with tab2:
#    st.warning('This feature is coming soon')
#    if 'dataset01' in st.session_state:
#        get_acts_perDs(st.session_state.dataset01)
with tab3: 
    st.warning('This feature is coming soon')







    ### Conditional Behaviour ###
    
    # Update params and Run
    #w_run = st.button("Run Analysis")
    #if w_run:     

    run_counter =  1
    
        
    #if w_run and (run_counter == 1):
        #tab5 = st.tabs("Results of Dataset1")
    #    with tab2:
    #        st.write(st.session_state.ds1)
        
with tab2:
    if 'genelist_ds1' in st.session_state:
        util.add_results(st.session_state['genelist_ds1'][0], st.session_state['genelist_ds1'][1], 'download_csv_ds1_sessionstate')


    













    #check = st.checkbox('check', value=False, key=None, help=None, on_change=None, args=None, kwargs=None, disabled=False, label_visibility="visible")
    #printcheck = st.checkbox('Print')
    #if printcheck:
    #            st.subheader('State')
    #            st.write("miau")
    #            st.text_input("input")
    #            st.subheader('Action')
    #            st.write("huhu")
    #            st.subheader('Reward')
    #            st.write("he")
    #if input == "x": 
    #    "Sooooo Cute!"





def process_sc():
    """ Single Cell Data """
    #subprocess.run(f'poetry run python analysis_params.py', shell=True)
    analysis_params = util.get_analysis_params()
    dc_dataset = sc_classes.Analysis.new_dataset(sc_classes.Baseanalysis, dcu.Decoupler)
    datasets = [
            ('01', 'scRNA', 'mouse', dc_dataset)
            ]
    params_path = path.abspath('./') # /Users/hanna/Documents/projects/SGUI/CTLA4/v00/analysis/")  TODO: make relative path
    analysis = sc_classes.Analysis(datasets, params_path)
    util.print_info(analysis, True)
    scl.analysis = analysis
    util.save_paths(analysis)
    data = 'Success'

def display_genelist_result(df_melt, fig):
    """ In second tab
    """
    util.add_results(df_melt, fig, 'download_csv_main')
    with tab2:
        util.add_results(df_melt, fig, 'download_csv_ds1')
    if 'genelist_ds1' not in st.session_state:
        st.session_state['genelist_ds1'] = [df_melt, fig]






# """
# ### Analysis ###
#     w_start = st.button("Start")
#     if w_start:
        

#         dc_dataset = sc_classes.Analysis.new_dataset(sc_classes.Baseanalysis, dcu.Decoupler) 
#         analysis = sc_classes.Analysis(datasets=[
#                     ('01', 'scRNA', 'mouse', dc_dataset)
#                     ], params_path = path.abspath("/Users/hanna/Documents/projects/FUNKI/CTLA4/v00/analysis/"))
#         util.print_info(analysis, True)

#         scl.analysis = analysis
#         util.save_paths(analysis)
#         #scl.get_acts()
#         "Activity calculation is done."
#         #scl.plot_umap()
#         #scl.get_mean_acts()
#         #scl.plot_mean_acts()

#         # Paths
#         figp = analysis.datasets[0].analysis_params['paths']['figpath']
#         #resp = analysis.datasets[0].analysis_params['paths']['resultpath']
#         prog = 'file:///' + figp + '/log/decoupler/progeny/300_mlm/umaps/umap.pdf'
#         doro = glob.glob(figp + '/log/decoupler/dorothea/abc_mlm/umaps/*')  
#         doro_mean =   glob.glob(figp + '/log/decoupler/dorothea/abc_mlm/mean_acts/*')  

#         #import pandas as pd
#         #df = pd.read_csv('file:///' + resp + '/log/decoupler/progeny/300_mlm/')
#         #print(df.to_string()) 

#         with st.sidebar:
#             util.show_dorothea_res(doro[0:10], colNumb=1)
            
#         util.show_dorothea_res(doro_mean, colNumb=1, counter = 0, mean_act=True)
#         util.gui_paths(analysis)




#         ### Progeny Results ###
#         progRes = st.button('Open Progeny results in new browser tab')
#         if progRes:
#             webbrowser.open_new_tab(prog)
#         st.success('Models are loaded!')

# """




# """
# UseCase1: gene list
#   case 1: default
#   case 2: test data
#   case 3: phospho
# UseCase2: matrix
#   case 1: default
#   case 2: test data
#   case 3: phospho
# UseCase3: h5ad
#   case 1: default
#   case 2: test data
# """


#add_nested_key(analysis_params, ['dataset_params', self.organism, self.seq_type])
    #organism = list(analysis_params['dataset_params'].keys())[0]
    #sc_funcs.dict_replace(analysis_params, w_organism, sc_funcs.getpath(organism))
    #analysis_params['dataset_params'][w_organism][w_omicstype] = analysis_params['dataset_params'].pop((list(analysis_params['dataset_params'])[0]))
   


#data_canada = px.data.gapminder().query("country == 'Canada'")    
#example_df.loc[example_df["column_name1"] condition, "column_name2"] = value
#scat.scatter(x, y ,s=s, label='p-value')
#plt.ylim(-1,1)
#plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5), labelspacing=3)

