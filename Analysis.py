import webbrowser, glob, os, subprocess              # json, yaml, base64, re, sys, numba, logging, random, dill, logging.config
from os import makedirs, path       
from copy import deepcopy
import streamlit as st, decoupler as dc, numpy as np, pandas as pd 
st.set_page_config(layout="wide")   ##import matplotlib.pyplot as plt
#from standard_workflows import default_analysis_params as ap
from standard_workflows import * 
from gui import pages, utils
from gui.utils import utilities as util
UiVal = util.UiVal

# st.success(f'The analysis will be run for **{w_organism}** with **{w_omicstype}** data in the format **{w_inputformat}**.')
# #st.session_state.genelist = "Please wait while your data is processing..."
#"""We use log10 for the plot because it has a better interpretability than log2.
#We can cache data of a function with the decorator: @st.cache_data . This way a function that does some heavy calculations doesn't recalculate when it gets the same input again but receives the data from the cache instead. 
#"""
# sparse matrix wrong datatype object instead of float in df coming from pd.read_csv
# row t, cols genes -> row t, cols pathways
# Ora background genes: we only need the number of background genes and not the actual genes but it's more intuitive for the people if they can just upload the genelist instead of counting the genes.
# Tabs: 'Prior Knowledge', 'COSMOS/CARNIVAL', 'Pseudobulk', 'Parameter Choice'
# #d.values.astype(float).index
# TODO: 
# - Same project and datasetname -> read in analysis params file if exists and execute analsis accordingly on 'submit'
# - Add new dataset
# - Choose already processed datasets that shall show their results in the tabs
# - keep plotting function very general to use it with varous input data (left col for table, right col for plot)

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



util.init_page()



def get_acts(w_matrix, group):
    """run decoupler

    Args:
        w_matrix (DataFrame): first column gene names, second column p-values or other statistics
        w_organism (_type_): _description_
        group (_type_): _description_
    """
    data = w_matrix
    ap = st.session_state.ap
    priorKnwldg = sc_funcs.deep_get(ap['dataset_params'], (sc_funcs.getpath(ap['dataset_params'], 'priorKnowledge')))['priorKnowledge'].keys()
    organism = list(ap['dataset_params'].keys())[0]
    for resource in priorKnwldg:
        param = sc_funcs.deep_get(ap, sc_funcs.getpath(ap, resource, search_value = False))
        if resource != 'collectri':
            net = eval(f'dc.get_{resource}(organism,' + str(list(param.values())[0][0]) + ')')
        else: 
            net = pd.read_csv(ap['proj_params']['paths']['data_root_path']+'/collectri.csv')
        d = data.transpose().tail(-1) 
        d.columns = data[data.columns[0]]
        d = d.astype(float)
        result = dc.decouple(d, net)
        method = 'ulm'
        result[f'{method}_pvals'].index = ['pvals']
        result_pvals = result[f'{method}_pvals'].transpose()
        result_pvals[resource] = result_pvals.index
        result_estimate = result[f'{method}_estimate'].transpose()
        result_estimate[resource] = result_estimate.index
        df = pd.merge(result_pvals, result_estimate)
        #df = df[['genes', 't', 'pvals']] # reorder cols
        df.index = df[resource]
        df = df.drop(resource, axis = 1)
        figpath = f'{resource}.png'
        dc.plot_barplot(result[f'{method}_estimate'], 't', top=25, vertical=False, return_fig = False, save = figpath)
        util.add_results(df, figpath, resource.capitalize(), resource)

def process_genelist(w_genelist, w_organism):
    """Process the gene names provided by the user. 
    Uses get_ora_df as it is the only method that can be used with a gene list as input instead of a matrix.
    """
    features = w_genelist
    data = pd.DataFrame({'group': ['undefined']*len(features), 'features': features})
    def calculate_prior(net, priortype):
        """ Get activities from decoupler

        Args:
            features (String): name of column with significant features
            net (DataFrame): prior knwoledge retrieved via decoupler 
            priortype (String): column name for result

        Returns:
            List: [DataFrame, Plotly, Title]
        """
        res = dc.get_ora_df(data, net, groupby='group', features='features')
        # melt to long format
        res = res.reset_index().melt(id_vars=['index'], var_name=priortype, value_name='p-value')
        # sort by the p-values and assign ranks
        res = res.sort_values(by='p-value')
        #res['rank'] = list(range(1, len(res.index) + 1))
        # Remove group, reset index, change column order
        res = res.drop('index', axis=1).reset_index(drop=True)
        res = res[[priortype, 'p-value']]
        #res = res[['rank', priortype, 'p-value']]
        import plotly.express as px
        import numpy as np
        # to avoid divide by zero error 
        #newp = np.where(res['p-value'] > 0.0000000001, res['p-value'] , -10)
        #res['log10(p-value)'] = np.log10(newp, out=newp, where=newp > 0)*-10
        res['log10(p-value)'] = np.log10(res['p-value'])*-10
        if len(res) >= 20: 
            max_rows = 20
        else: 
            max_rows = len(res)
        res_plot = res.iloc[0:max_rows, :]
        if(res_plot['log10(p-value)'].value_counts().max() == 20):
            print(res_plot)
            fig = 'false'
        else: 
            fig = px.bar(res_plot, x=priortype, y='log10(p-value)')
        return [res.iloc[:, 0:2], fig, priortype.capitalize()] # result, figure, title

    prog = calculate_prior(dc.get_progeny(organism = w_organism), 'pathway')
    collectri = calculate_prior(pd.read_csv(st.session_state.ap['proj_params']['paths']['data_root_path']+'/collectri.csv'), 'transcriptionFactor')
    util.add_results(prog[0], prog[1], prog[2],'res_prog')
    util.add_results(collectri[0], collectri[1], collectri[2],'res_collectri')


def get_testdata(w_inputformat, analysis_params):
    """Read and process Testdata"""
    data = ''
    match w_inputformat:
        case UiVal.GENES:
            data = ['KIAA0907', 'KDM5A', 'CDC25A', 'EGR1', 'GADD45B', 'RELB', 'TERF2IP', 'SMNDC1', 'TICAM1', 'NFKB2', 'RGS2', 'NCOA3', 'ICAM1', 'TEX10', 'CNOT4', 'ARID4B', 'CLPX', 'CHIC2', 'CXCL2', 'FBXO11', 'MTF2', 'CDK2', 'DNTTIP2', 'GADD45A', 'GOLT1B', 'POLR2K', 'NFKBIE', 'GABPB1', 'ECD', 'PHKG2', 'RAD9A', 'NET1', 'KIAA0753', 'EZH2', 'NRAS', 'ATP6V0B', 'CDK7', 'CCNH', 'SENP6', 'TIPARP', 'FOS', 'ARPP19', 'TFAP2A', 'KDM5B', 'NPC1', 'TP53BP2', 'NUSAP1', 'SCCPDH', 'KIF20A', 'FZD7', 'USP22', 'PIP4K2B', 'CRYZ', 'GNB5', 'EIF4EBP1', 'PHGDH', 'RRAGA', 'SLC25A46', 'RPA1', 'HADH', 'DAG1', 'RPIA', 'P4HA2', 'MACF1', 'TMEM97', 'MPZL1', 'PSMG1', 'PLK1', 'SLC37A4', 'GLRX', 'CBR3', 'PRSS23', 'NUDCD3', 'CDC20', 'KIAA0528', 'NIPSNAP1', 'TRAM2', 'STUB1', 'DERA', 'MTHFD2', 'BLVRA', 'IARS2', 'LIPA', 'PGM1', 'CNDP2', 'BNIP3', 'CTSL1', 'CDC25B', 'HSPA8', 'EPRS', 'PAX8', 'SACM1L', 'HOXA5', 'TLE1', 'PYGL', 'TUBB6', 'LOXL1']
            process_genelist(data, w_organism)
            #display_genelist_result(result[0], result[1])
        case UiVal.MATRIX:
            data = pd.read_csv(st.session_state.ap['proj_params']['paths']['data_root_path']+'/differential_stats.csv')
            get_acts(data, data.columns[1])
        case UiVal.H5AD: 
            process_sc()
    return data


def get_data(w_inputformat, analysis_params):
    match w_inputformat:
        case UiVal.GENES:
            w_genelist = st.text_area("Please paste your list of comma separated gene names (i.e. DEGs) here:", key= "genelist", placeholder= "Gene1, Gene2") #on_change=lambda x:send_genelist(x), args=(st.session_state["genelist"]))
            w_analyse_genes = st.button('Analyse Genes')
            if w_analyse_genes:
                process_genelist(w_genelist.split(', '), w_organism)
                #display_genelist_result(result[0], result[1])
        case UiVal.MATRIX:
            st.caption('Please provide a csv file with the gene names in the first column of the table.')
            uploaded_files = st.file_uploader("Choose a CSV file", accept_multiple_files=True, type = ".csv") # TODO: allow tsv
            if len(uploaded_files) >= 1:
                cols = st.columns(len(uploaded_files)) 
                for i in range(0, len(uploaded_files)):
                    file = uploaded_files[i]
                    data = pd.read_csv(file)
                    with cols[i]:
                        st.write(file.name, data)
                    if(data.shape[1] == 1):
                        process_genelist(data[1:].to_csv(header=None, index=False).strip('\n').split('\n'), w_organism)
                    else:
                        get_acts(data, data.columns[1])
        case UiVal.H5AD:
            sc.success('Congrats, your data unlocks the project management feature! This is an additional service that automatically downloads all results in a reproducible way.')
            #"/Users/hanna/Documents/projects/SGUI/CTLA4/v00/analysis/mouse/scRNA/01/data/01.h5ad"
            projpath = st.text_input('basepath')
            projname = st.text_input('project name')
            input_path = st.text_input('input data path')
            ok = st.button('OK')
            if w_testdata:
                projpath = './'
                projname = w_projname
                input_path = './example_inputs/'
            if ok & (projname != None) & (projpath != None):
                tbdatadir = path.join(projpath, projname, 'v00/analysis/01/data/')
                tbdatapath = path.join(tbdatadir, '01.h5ad')
                
                if not path.exists(tbdatadir):
                    os.makedirs(tbdatadir)
                import scanpy as sc 
                sc.write(tbdatapath, sc.read(input_path, cache = True))

            projpath = path.join(projpath, projname)
            scriptpath = path.join(projpath, 'scripts/python/')
            if not path.exists(scriptpath):
                os.makedirs(scriptpath)
            subprocess.run(f'cp /Users/hanna/Documents/projects/SGUI/analysis_params.py {scriptpath}/analysis_params.py', shell=True)


st.header("Transcription Factor And Pathway Analysis")
st.markdown("**Analysis of bulk RNA sequencing data with the tool Decoupler**")


tab1, tab2, tab3, tab4 = st.tabs(['Analysis', 'Results 1st dataset', 'Results 2nd dataset', 'Parameter Choices'])

# TODO: if not 'Run', don't change any params or results.
# TODO: Extend to multiple datasets at once? Not really necessary. Maybe one tab per ds results?

w_datasetname_default = '01'
with tab1:
    ### UI ###

    # Dataset Specific Params
    st.caption('Describe your data')
    paramcol1, paramcol2 = st.columns(2)
    with paramcol1:
        w_organism    = st.selectbox('Organism', (UiVal.HUMAN, UiVal.MOUSE))
        w_omicstype   = st.selectbox('Omics Type', (UiVal.BULKRNA, UiVal.PHOSPHO, UiVal.SCRNA,)) 
        w_inputformat = st.selectbox('Input format', (UiVal.GENES, UiVal.MATRIX, UiVal.H5AD))
        # testdata
        w_testdata    = st.checkbox('Use test data')
    st.session_state.ap = util.get_analysis_params(w_organism, w_omicstype)
    if not w_testdata:
        get_data(w_inputformat, st.session_state.ap)
    else:
        data = get_testdata(w_inputformat, st.session_state.ap)
    with paramcol2:
        # Project Specific Params
        w_show_all_opts = st.checkbox('Show all options')
        if(w_show_all_opts): 
            st.caption('Advanced Settings')
            w_projname_default = st.session_state.ap['proj_params']['proj_id']
            w_projname    = st.text_input(UiVal.PROJ, placeholder = w_projname_default)
            if w_projname != '':
                util.change_param('proj_id', w_projname)
            else:
                w_projname = st.session_state.ap['proj_params']['proj_id']
            #datasetnames   = list(st.session_state.ap['dataset_params'][w_organism][w_omicstype].keys())
            #datasetnames.pop() # drop the prior knowledge element
            
            w_datasetname = st.text_input('Dataset Name', placeholder = w_datasetname_default)#datasetnames)
            if w_datasetname == '':
                w_datasetname = w_datasetname_default
            #TODO: w_resultpath  = st.text_input("Result Path (if empty, results won't be saved on your computer)", './')    
            #TODO: topn = st.text_input('top n', '')
            #TODO: st.write('You chose the following parameters for topN: ', list(topn.split(',')))
            st.success(f'The results will be saved in **./{w_projname}/{w_datasetname}/**.')

            sc_funcs.dict_delete_key(st.session_state.ap, sc_funcs.getpath(st.session_state.ap, w_datasetname_default))
            st.session_state.ap['dataset_params'][w_organism][w_omicstype][w_datasetname] = {}
    
    with tab4:
        st.write('### Chosen Analysis Parameters')
        st.write(st.session_state.ap)
        ap = st.session_state.ap
        st.write('### Merged')
        if 'w_datasetname' not in locals():
            w_datasetname = w_datasetname_default
        st.write(sc_funcs.merge_dicts(st.session_state.ap['proj_params'], st.session_state.ap['dataset_params'][w_organism][w_omicstype][w_datasetname]) )
        st.write('### Session State')
        st.write(st.session_state)
        st.write('### Paths')
        #st.markdown(sc_funcs.print_paths(st.session_state.ap['proj_params']['paths']))

    

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

