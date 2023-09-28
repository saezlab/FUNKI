import streamlit as st
import os, base64, re, sys
from os import makedirs, path
#from copy import deepcopy
#from IPython.display import display, Markdown   #, Latex # to display Markdown in code chunk
# json, yaml, numba, logging, random, dill, logging.config
from standard_workflows import *  
sys.path.append('../')
from streamlit_extras.switch_page_button import switch_page
from st_pages import Page, show_pages, _hide_pages
from gui.utils import utilities as util
import pandas as pd
import decoupler as dc
#import openpyxl
import logging
import plotly.express as px
import numpy as np
import kaleido
import subprocess
import glob
#########################
# TODO: allow tsv input #
#########################

UiVal = util.UiVal
#@st.cache_data(experimental_allow_widgets=True)
def get_acts(dataset):
    """run decoupler

    Args:
        w_matrix (DataFrame): first column gene names, second column p-values or other statistics
        w_organism (_type_): _description_
        group (_type_): _description_
    """
    data = dataset['data']
    ap = st.session_state.ap

    priorKnwldg = sc_funcs.deep_get(ap['dataset_params'], (sc_funcs.getpath(ap['dataset_params'], 'priorKnowledge')))['priorKnowledge'].keys()
    organism = list(ap['dataset_params'].keys())[0]
    for resource_type in priorKnwldg:
        resources = sc_funcs.deep_get(ap['dataset_params'], (sc_funcs.getpath(ap['dataset_params'], resource_type)))[resource_type].keys()
        datarootpath = st.session_state.ap['proj_params']['paths']['data_root_path']

        #subprocess.call([f'{datarootpath}getfrom_omnipath.R {resource} {organism}'], shell = True)
        for resource in resources:
            param = sc_funcs.deep_get(ap, sc_funcs.getpath(ap, resource, search_value = False))
            
            net = pd.DataFrame()
            if resource == 'progeny' and organism == 'mouse':
                    param_val = str(list(param.values())[0][0])
                    priorKnwldg_avail_path = r'./data/priorKnowledge/*' + re.escape(resource) + '*' + re.escape(organism) + '*' + re.escape(param_val) + '*.csv'
                    priorKnwldg_avail = glob.glob(priorKnwldg_avail_path)
                    if priorKnwldg_avail == [] :
                        priorKnwldg_avail_path = r'./data/priorKnowledge/*' + re.escape(resource) + '*' + re.escape(organism) + '*.csv'
                    priorKnwldg_avail = glob.glob(priorKnwldg_avail_path)[0]
                    #priorKnwldg_avail = [os.path.basename(x) for x in glob.glob(priorKnwldg_avail_path)]
                    #st.write(f'The following file is used as priorKnowledge resource:  \n\n`{priorKnwldg_avail}`')
                    net = pd.read_csv(priorKnwldg_avail, index_col = False)
            elif resource == 'ksn':
                if organism == 'human':
                    net = dc.get_ksn_omnipath()
                else: 
                    st.warning('This priorKnowledgeResource is only available for human data so far.')
                    return
            else:
                # st.write(f'**PriorKnowledge Resource {resource}**')
                print(f'{resource}, {organism}, {str(list(param.values())[0][0])}')
                net = eval(f'dc.get_{resource}(organism,' + str(list(param.values())[0][0]) + ')')
            
            #st.write(net)
            method = dataset['method']
            figpath = f'{resource}.png'
            title = f'{resource_type.capitalize()} retrieved from {resource.capitalize()}'
            if method == 'ora_df':
                d = data
                if type(d) == pd.DataFrame:
                    d = d[d.columns[0]].tolist()

                result = dc.get_ora_df(d, net)
                def plot_pval(data, idCol, pvalCol, figpath, max_rows = 20):
                    """barplot showing -log10 pvals but not more than 20 results

                    Args:
                        data (DataFrame): full dataframe
                        idCol (str): x axis
                        pvalCol (str): y axis, gets transformed with -log10
                        figpath (str): path where the resulting figure is saved

                    Returns:
                        _type_: _description_
                    """
                    data = data[[idCol, pvalCol]].sort_values(by=pvalCol)
                    # assign ranks
                    data.columns = [resource, pvalCol]
                    data[f'-log10({pvalCol})'] = np.log10(data[pvalCol])*-1
                    if len(result) < max_rows: 
                        max_rows = len(data)
                    data_plot = data.iloc[0:max_rows, :]
                    
                    if(data_plot[f'-log10({pvalCol})'].value_counts().max() == 20):
                        print(data_plot)
                        fig = 'false'
                    else: 
                        fig = px.bar(data_plot, x=resource, y='-log10(p-value)')
                        fig.write_image(figpath)
                    return data, fig
                df, fig = plot_pval(result, 'Term', 'p-value', figpath)                
            else:
                try:
                    d = data.transpose().tail(-1)  # drop first row
                    d.columns = data[data.columns[0]] # reset column names
                    d = d.astype(float)

                    result = dc.decouple(d, net, methods=[method])
                    print(result)
                    def decouple_result_totable(result):
                        result[f'{method}_pvals'].index = ['pvals']
                        
                        result_pvals = result[f'{method}_pvals'].transpose()
                        result_pvals[resource] = result_pvals.index
                        
                        result_estimate = result[f'{method}_estimate'].transpose()
                        result_estimate[resource] = result_estimate.index
                        result_estimate.columns = ['t', resource] # t is estimated activity

                        df = pd.merge(result_pvals, result_estimate)
                        df.index = df[resource]
                        df = df.drop(resource, axis = 1)
                        return df
                    df = decouple_result_totable(result)
                    
                    dc.plot_barplot(result[f'{method}_estimate'], result[f'{method}_estimate'].index[0], top=25, vertical=False, return_fig = False, save = figpath)
                    
                except ValueError as ve:
                    try:
                        result = dc.decouple(d, net, methods=[method], min_n=2)
                        st.warning("Warning: All sources with at least two targets were taken into consideration. The default would be to have at least five targets per source. To go with the default you would need to add more genes.")
                    except ValueError as ve:
                        st.error("Error: There aren't any sources with the minimum of five targets. Please provide more genes.")
            util.add_results(df, figpath, title, f"{title}{dataset['datasetid']}")
            


def get_data(w_inputformat)->list[dict] :
    datasets = list()
    def increment_analysiscount():
        st.session_state.analysiscount += 1
    match w_inputformat:
        case UiVal.GENES | UiVal.KINASES:
            w_elementlist = st.text_area("Please paste your list of comma separated element names (i.e. DEGs,...) here:", key= "elementlist", placeholder= "element1, element2") #on_change=lambda x:send_genelist(x), args=(st.session_state["genelist"]))
            data = w_elementlist.split(', ')
            #data = data.split(',')
            #data = data.split('\t')
            w_analyse_elements = st.button('Analyse', key=f'analyse_button', on_click=increment_analysiscount)
            if w_analyse_elements:
                datasets.append({'datasetname': '01', 'data': data, 'method': 'ora_df', 'datasetid': 1})
            st.write('You can find the results in the "dataset1" tab above.')
        case UiVal.MATRIX: # TODO add warning when NA included
            st.caption('Please provide a csv, tsv or xlsx file.')
            uploaded_files = st.file_uploader("Choose a file", accept_multiple_files=True, type = ['.xlsx', '.csv']) 
            
            # read and display data per upload and return 'datasets'
            if len(uploaded_files) >= 1:
                def read_and_display_files(uploaded_files)->list[dict]:
                    datasets = list()
                    data = ''
                    method = ''
                    st.write('Please provide the names of your entity column and the statistics column if available. Leave the stats column name empty if your data has no statistics or consists of significant elements only. By default the fields are filled with the names of the first and second column of the uploaded data table.')
                    
                    
                    w_analyse_elements = st.button('Analyse', key=f'analyse_button', on_click=increment_analysiscount)

                    
                    filecols = st.columns(len(uploaded_files)) 
                    for i in range(0, len(uploaded_files)):
                        # read data
                        file = uploaded_files[i]
                        filename = os.path.splitext(file.name)[0].lower()
                        data = util.read_file(file)

                        with filecols[i]:
                            # display data
                            w_choose_id_col, w_choose_stats_col = st.columns(2)
                            with w_choose_id_col:
                                id_colname = st.selectbox('ID column', data.columns, 0, key=f'idcolname_{i}')
                                
                            with w_choose_stats_col:
                                if len(data.columns) >= 2:
                                    cols = list(data.columns)
                                    cols.append(None)
                                    stats_colname = st.selectbox('Statistics column', cols, 1, key=f'statscolname_{i}')
                                else:
                                    stats_colname = st.selectbox('Statistics column', list() ,0, key=f'statscolname_{i}')

                            st.write(f'Dataset{i}: {filename}', data)

                        
                            # return relevant information per dataset
                            if(w_analyse_elements):    
                                # clean data
                                if(stats_colname == None):
                                    data = data[[id_colname]]
                                    method = 'ora_df'
                                else:
                                    data = data[[id_colname, stats_colname]]
                                    method = 'ulm'
                                data[id_colname] = data[id_colname].apply(str.upper)
                                datasets.append({'datasetname': filename, 'data': data, 'method': method, 'datasetid': i})
                                st.write('You can find the results in the "dataset1" tab above.')
                    return datasets
                            
                datasets = read_and_display_files(uploaded_files)
    return datasets




def get_testdata(w_inputformat, w_omicstype, datarootpath):
    """Read and process Testdata"""
    data = ''
    datasets = list()
    match w_inputformat:
        case UiVal.KINASES:
            data = ['Tex2_S265', 'Abl2_S822', 'Abl2_S819', 'Prkd3_S213', 'Tfeb_S108', 'Prkd3_S216', 'Foxk2_S389', 'Mtor_T2473', 'Peak1_S1203', 'Peak1_S1206', 'Lsr_S591', 'Eps15l1_S253', 'Lsr_S588', 'Dennd6a_S16', 'Unc5b_S531', 'Nfia_S310', 'Mapk6_S683', 'Mapk6_S687', 'Unc5b_S528', 'Tns3_Y773', 'Tfeb_S121', 'Mybbp1a_T1260', 'Washc2_S1048', 'Washc2_S1049', 'Tfeb_S113', 'Clip1_S194', 'Iws1_S667', 'Larp1_S743', 'Znf608_S626', 'Trmt10c_S79', 'Gpatch8_S733', 'Gpatch8_S735', 'Fnbp1_S299', 'Nfia_S303', 'Ankrd17_S1692', 'Rhbdf1_T180', 'Pnisr_S321', 'Trmt10c_S85', 'Micall1_S494', 'Micall1_S496', 'Ankrd17_S1696', 'Top2b_S1448', 'Rabep1_S410', 'Znf608_T635', 'Tnks1bp1_S889', 'Mrpl11_S160', 'Nufip2_T221', 'Morc2a_S737', 'Rgl3_S568', 'Rgl3_S572', 'Ercc5_S572', 'Lsm11_S18', 'Arhgap12_S211', 'Rbm25_S678', 'Samd4b_S585', 'Sipa1l1_S211', 'Supt5h_T661', 'Srrm2_S2691', 'Srrm2_T2689', 'Plec_S4396', 'Iws1_S666', 'Pnn_S698', 'Pnn_S700', 'Lmna_T19', 'Phldb2_S71']
            elementtype = 'kinases'
            elements = "Tex2_S265, Abl2_S822, Abl2_S819, Prkd3_S213, Tfeb_S108, Prkd3_S216, Foxk2_S389, Mtor_T2473, Peak1_S1203, Peak1_S1206, Lsr_S591, Eps15l1_S253, Lsr_S588, Dennd6a_S16, Unc5b_S531, Nfia_S310, Mapk6_S683, Mapk6_S687, Unc5b_S528, Tns3_Y773, Tfeb_S121, Mybbp1a_T1260, Washc2_S1048, Washc2_S1049, Tfeb_S113, Clip1_S194, Iws1_S667, Larp1_S743, Znf608_S626, Trmt10c_S79, Gpatch8_S733, Gpatch8_S735, Fnbp1_S299, Nfia_S303, Ankrd17_S1692, Rhbdf1_T180, Pnisr_S321, Trmt10c_S85, Micall1_S494, Micall1_S496, Ankrd17_S1696, Top2b_S1448, Rabep1_S410, Znf608_T635, Tnks1bp1_S889, Mrpl11_S160, Nufip2_T221, Morc2a_S737, Rgl3_S568, Rgl3_S572, Ercc5_S572, Lsm11_S18, Arhgap12_S211, Rbm25_S678, Samd4b_S585, Sipa1l1_S211, Supt5h_T661, Srrm2_S2691, Srrm2_T2689, Plec_S4396, Iws1_S666, Pnn_S698, Pnn_S700, Lmna_T19, Phldb2_S71"
            #data = pd.DataFrame({'features': data})
            st.markdown(f"**The following {elementtype} are used:**")
            st.write(elements)
            datasets = {'datasetname': 'genesymbol_phosphosites', 'data': data, 'method': 'ora_df', 'datasetid': 0}
        case UiVal.GENES:
            data = ['KIAA0907', 'KDM5A', 'CDC25A', 'EGR1', 'GADD45B', 'RELB', 'TERF2IP', 'SMNDC1', 'TICAM1', 'NFKB2', 'RGS2', 'NCOA3', 'ICAM1', 'TEX10', 'CNOT4', 'ARID4B', 'CLPX', 'CHIC2', 'CXCL2', 'FBXO11', 'MTF2', 'CDK2', 'DNTTIP2', 'GADD45A', 'GOLT1B', 'POLR2K', 'NFKBIE', 'GABPB1', 'ECD', 'PHKG2', 'RAD9A', 'NET1', 'KIAA0753', 'EZH2', 'NRAS', 'ATP6V0B', 'CDK7', 'CCNH', 'SENP6', 'TIPARP', 'FOS', 'ARPP19', 'TFAP2A', 'KDM5B', 'NPC1', 'TP53BP2', 'NUSAP1', 'SCCPDH', 'KIF20A', 'FZD7', 'USP22', 'PIP4K2B', 'CRYZ', 'GNB5', 'EIF4EBP1', 'PHGDH', 'RRAGA', 'SLC25A46', 'RPA1', 'HADH', 'DAG1', 'RPIA', 'P4HA2', 'MACF1', 'TMEM97', 'MPZL1', 'PSMG1', 'PLK1', 'SLC37A4', 'GLRX', 'CBR3', 'PRSS23', 'NUDCD3', 'CDC20', 'KIAA0528', 'NIPSNAP1', 'TRAM2', 'STUB1', 'DERA', 'MTHFD2', 'BLVRA', 'IARS2', 'LIPA', 'PGM1', 'CNDP2', 'BNIP3', 'CTSL1', 'CDC25B', 'HSPA8', 'EPRS', 'PAX8', 'SACM1L', 'HOXA5', 'TLE1', 'PYGL', 'TUBB6', 'LOXL1']
            elementtype = 'genes'
            elements = "KIAA0907, KDM5A, CDC25A, EGR1, GADD45B, RELB, TERF2IP, SMNDC1, TICAM1, NFKB2, RGS2, NCOA3, ICAM1, TEX10, CNOT4, ARID4B, CLPX, CHIC2, CXCL2, FBXO11, MTF2, CDK2, DNTTIP2, GADD45A, GOLT1B, POLR2K, NFKBIE, GABPB1, ECD, PHKG2, RAD9A, NET1, KIAA0753, EZH2, NRAS, ATP6V0B, CDK7, CCNH, SENP6, TIPARP, FOS, ARPP19, TFAP2A, KDM5B, NPC1, TP53BP2, NUSAP1, SCCPDH, KIF20A, FZD7, USP22, PIP4K2B, CRYZ, GNB5, EIF4EBP1, PHGDH, RRAGA, SLC25A46, RPA1, HADH, DAG1, RPIA, P4HA2, MACF1, TMEM97, MPZL1, PSMG1, PLK1, SLC37A4, GLRX, CBR3, PRSS23, NUDCD3, CDC20, KIAA0528, NIPSNAP1, TRAM2, STUB1, DERA, MTHFD2, BLVRA, IARS2, LIPA, PGM1, CNDP2, BNIP3, CTSL1, CDC25B, HSPA8, EPRS, PAX8, SACM1L, HOXA5, TLE1, PYGL, TUBB6, LOXL1"
            #data = pd.DataFrame({'features': data})
            st.markdown(f"**The following {elementtype} are used:**")
            st.write(elements)
            datasets = {'datasetname': 'genelist', 'data': data, 'method': 'ora_df', 'datasetid': 0}
        case UiVal.MATRIX:
            match w_omicstype:
                case UiVal.BULKRNA:
                    data = pd.read_csv(datarootpath+'/differential_stats.csv')
                    datasets= {'datasetname': 'differential_stats', 'data': data, 'method': 'ulm', 'datasetid': 0}
                case UiVal.PHOSPHO:
                    data = pd.read_excel(datarootpath +'/dePhosRes.xlsx')
                    data = data.loc[:,['site','t_statistic']]
                    data['site'] = data['site'].apply(lambda x: x.upper())
                    data['site'] = data['site'].apply(str.upper)
                    datasets= {'datasetname': 'smartphos', 'data': data, 'method': 'ulm', 'datasetid': 0}
        case UiVal.H5AD: 
            process_sc()
        
    return datasets


