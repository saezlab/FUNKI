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
#########################
# TODO: allow tsv input #
#########################

UiVal = util.UiVal
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
        resources = resources = sc_funcs.deep_get(ap['dataset_params'], (sc_funcs.getpath(ap['dataset_params'], resource_type)))[resource_type].keys()
        for resource in resources:
            param = sc_funcs.deep_get(ap, sc_funcs.getpath(ap, resource, search_value = False))
            if resource == 'collectri':
                datarootpath = st.session_state.ap['proj_params']['paths']['data_root_path']
                net = pd.read_csv(f"{datarootpath}/collectri.csv")
            elif resource != 'ksn':
                print(f'{resource}, {organism}, {str(list(param.values())[0][0])}')
                net = eval(f'dc.get_{resource}(organism,' + str(list(param.values())[0][0]) + ')')
            else: 
                #net = pd.read_csv(ap['proj_params']['paths']['data_root_path']+'/collectri.csv')
                net = dc.get_ksn_omnipath()
                print(net)
            d = data.transpose().tail(-1) 
            d.columns = data[data.columns[0]]
            d = d.astype(float)
            try:
                print(net)
                print(d)
                d.to_csv("./input.csv")
                method = dataset['method']
                if method == 'ora_df':
                    st.write('this function is work in progress')
                    return
                else:
                    result = dc.decouple(d, net, methods=[method])
                #st.write(result)
                print(result)
                #st.write(result[f'{method}_pvals'])
                result[f'{method}_pvals'].index = ['pvals']
                result_pvals = result[f'{method}_pvals'].transpose()
                result_pvals[resource] = result_pvals.index
                result_estimate = result[f'{method}_estimate'].transpose()
                result_estimate[resource] = result_estimate.index
                result_estimate.columns = ['t', resource]
                df = pd.merge(result_pvals, result_estimate)
                #df.columns = [resource, 'pvals', 't']
                #df = df[['genes', 't', 'pvals']] # reorder cols
                df.index = df[resource]
                df = df.drop(resource, axis = 1)
                #st.write(df)
                figpath = f'{resource}.png'
                dc.plot_barplot(result[f'{method}_estimate'], result[f'{method}_estimate'].index[0], top=25, vertical=False, return_fig = False, save = figpath)
                title = f'{resource_type.capitalize()} retrieved from {resource.capitalize()}'
                util.add_results(df, figpath, title, f"{title}{dataset['datasetid']}")
            except ValueError as ve:
                try:
                    result = dc.decouple(d, net, min_n=2)
                    st.warning("Warning: All sources with at least two targets were taken into consideration. The default would be to have at least five targets per source. To go with the default you would need to add more genes.")
                except ValueError as ve:
                    st.error("Error: There aren't any sources with the minimum of five targets. Please provide more genes.")
    

def process_elementlist(w_genelist, w_organism):
    """Process the gene names provided by the user. 
    Uses get_ora_df as it is the only method that can be used with a gene list as input instead of a matrix.
    """
    features = w_genelist
    data = pd.DataFrame({'features': features})
    def calculate_prior(net, priortype):
        """ Get activities from decoupler

        Args:
            features (String): name of column with significant features
            net (DataFrame): prior knwoledge retrieved via decoupler 
            priortype (String): column name for result

        Returns:
            List: [DataFrame, Plotly, Title]
        """
        logging.info(data)
        print(data)
        print(net)
        res = dc.get_ora_df(data, net)
        # melt to long format
        res = res.reset_index().melt(id_vars=['index'], var_name=priortype, value_name='p-value')
        # sort by the p-values and assign ranks
        res = res.sort_values(by='p-value')
        #res['rank'] = list(range(1, len(res.index) + 1))
        # reset index, change column order
        res = res[[priortype, 'p-value']]
        #res = res[['rank', priortype, 'p-value']]
        import plotly.express as px
        import numpy as np
        # to avoid divide by zero error 
        #newp = np.where(res['p-value'] > 0.0000000001, res['p-value'] , -10)
        #res['log10(p-value)'] = np.log10(newp, out=newp, where=newp > 0)*-10
        res['-log10(p-value)'] = np.log10(res['p-value'])*-1
        if len(res) >= 20: 
            max_rows = 20
        else: 
            max_rows = len(res)
        res_plot = res.iloc[0:max_rows, :]
        if(res_plot['-log10(p-value)'].value_counts().max() == 20):
            print(res_plot)
            fig = 'false'
        else: 
            fig = px.bar(res_plot, x=priortype, y='-log10(p-value)')
        return [res.iloc[:, 0:2], fig, priortype.capitalize()] # result, figure, title
    try:
        prog = calculate_prior(dc.get_progeny(organism = w_organism), 'pathway')
        collectri = calculate_prior(pd.read_csv(st.session_state.ap['proj_params']['paths']['data_root_path']+'/collectri.csv'), 'transcriptionFactor')
        util.add_results(prog[0], prog[1], prog[2],'res_prog')
        util.add_results(collectri[0], collectri[1], collectri[2],'res_collectri')
    except ValueError as ve:
        try:
            prog = calculate_prior(dc.get_progeny(organism = w_organism, min_n = 2), 'pathway')
            collectri = calculate_prior(pd.read_csv(st.session_state.ap['proj_params']['paths']['data_root_path']+'/collectri.csv', organism = w_organism, min_n = 2), 'transcriptionFactor')
            util.add_results(prog[0], prog[1], prog[2],'res_prog')
            util.add_results(collectri[0], collectri[1], collectri[2],'res_collectri')
            st.warning("Warning: All sources with at least two targets were taken into consideration. The default would be to have at least five targets per source. To go with the default you would need to add more genes.")
        except ValueError as ve:
            st.error("Error: There aren't any sources with the minimum of five targets. Please provide more genes.")

def get_listdata(organism):
    """Provide a text field for string list input"""
    w_elementlist = st.text_area("Please paste your list of comma separated element names (i.e. DEGs,...) here:", key= "elementlist", placeholder= "element1, element2") #on_change=lambda x:send_genelist(x), args=(st.session_state["genelist"]))
    
    if w_analyse_elements:            
        process_elementlist(w_elementlist.split(', '), organism)
        #display_genelist_result(result[0], result[1])

def get_matrixdata(organims):
    st.caption('Please provide a csv or xlsx file.')
    uploaded_files = st.file_uploader("Choose a file", accept_multiple_files=True, type = ['.xlsx', '.csv']) 
            

# def get_data(w_inputformat, analysis_params):
#     data = ''
#     organism = list(analysis_params['dataset_params'].keys())[0]
#     match w_inputformat:
#         case UiVal.GENES:
#             w_elementlist = st.text_area(f"Please paste your list comma separated {UiVal.GENES} (i.e. DEGs) here:", key= "genelist", placeholder= "Gene1, Gene2") #on_change=lambda x:send_genelist(x), args=(st.session_state["genelist"]))
#             w_analyse_elements = st.button('Analyse')
#             if w_analyse_elements:            
#                 process_elementlist(w_elementlist.split(', '), organism)
#                 #display_genelist_result(result[0], result[1])
#         case UiVal.KINASES:
#                 ""
#         case UiVal.MATRIX: 
#             st.caption('Please provide a csv or xlsx file.')
#             uploaded_files = st.file_uploader("Choose a file", accept_multiple_files=True, type = ['.xlsx', '.csv']) 
#             if len(uploaded_files) >= 1:
#                 cols = st.columns(len(uploaded_files)) 
#                 for i in range(0, len(uploaded_files)):
#                     file = uploaded_files[i]
#                     ext = os.path.splitext(file.name)[-1].lower()
#                     match ext:
#                         case UiVal.CSV:
#                             data = pd.read_csv(file)
#                         case UiVal.EXCEL:
#                             data = pd.read_excel(file)

#                     with cols[i]:
#                         st.write(file.name, data)

#                     if(data.shape[1] == 1):
#                         process_genelist(data[1:].to_csv(header=None, index=False).strip('\n').split('\n'), w_organism)
#                     else:
#                         get_acts(data.iloc[:, 0:2], data.columns[1])
#         case UiVal.H5AD:
#             st.success('Congrats, your data unlocks the project management feature! This is an additional service that automatically downloads all results in a reproducible way.')
#             #"/Users/hanna/Documents/projects/SGUI/CTLA4/v00/analysis/mouse/scRNA/01/data/01.h5ad"
#             projpath = st.text_input('basepath')
#             projname = st.text_input('project name')
#             input_path = st.text_input('input data path')
#             ok = st.button('OK')
#             if w_testdata:
#                 projpath = './'
#                 projname = w_projname
#                 input_path = './example_inputs/'
#             if ok & (projname != None) & (projpath != None):
#                 tbdatadir = path.join(projpath, projname, 'v00/analysis/01/data/')
#                 tbdatapath = path.join(tbdatadir, '01.h5ad')
                
#                 if not path.exists(tbdatadir):
#                     os.makedirs(tbdatadir)
#                 import scanpy as sc 
#                 sc.write(tbdatapath, sc.read(input_path, cache = True))

#             projpath = path.join(projpath, projname)
#             scriptpath = path.join(projpath, 'scripts/python/')
#             if not path.exists(scriptpath):
#                 os.makedirs(scriptpath)
#             subprocess.run(f'cp /Users/hanna/Documents/projects/SGUI/analysis_params.py {scriptpath}/analysis_params.py', shell=True)



def get_data(w_inputformat):
    #filename = 'dataset01'
    data = ''
    method = ''
    i = 1
    datasets = list()
    match w_inputformat:
        case UiVal.GENES | UiVal.KINASES:
            w_elementlist = st.text_area("Please paste your list of comma separated element names (i.e. DEGs,...) here:", key= "elementlist", placeholder= "element1, element2") #on_change=lambda x:send_genelist(x), args=(st.session_state["genelist"]))
            data = w_elementlist.split(', ')
        case UiVal.MATRIX: # TODO add warning when NA included
            st.caption('Please provide a csv or xlsx file.')
            uploaded_files = st.file_uploader("Choose a file", accept_multiple_files=True, type = ['.xlsx', '.csv']) 
            data = ""
            # read and display data
            if len(uploaded_files) >= 1:
                st.write('Please provide the names of your entity column and the statistics column if available. Leave the stats column name empty if your data has no statistics or consists of significant elements only. By default the fields are filled with the names of the first and second column of the uploaded data table.')
                w_analyse_elements = st.button('Analyse', key=f'analyse_button')
                cols = st.columns(len(uploaded_files)) 
                for i in range(0, len(uploaded_files)):
                    file = uploaded_files[i]
                    ext = os.path.splitext(file.name)[-1].lower()
                    filename = os.path.splitext(file.name)[0].lower()
                    # read data
                    match ext:
                        case UiVal.CSV:
                            data = pd.read_csv(file)
                        case UiVal.EXCEL:
                            data = pd.read_excel(file)
                    # display data
                    with cols[i]:
                        columncol1, columncol2 = st.columns(2)
                        with columncol1:
                            id_guess = data.columns[0]   #data.iloc[:, 0:1] # guess that the first colum is the id column
                            id_colname = st.text_input('ID column', id_guess, key=f'idcolname_{i}')
                        with columncol2:
                            stats_guess = ''
                            if len(data.columns) >= 2:
                                stats_guess = data.columns[1] 
                            stats_colname = st.text_input('Statistics column', stats_guess, key=f'statscolname_{i}')
                       
                        st.write(f'Dataset{i}: {filename}', data)

                        if(w_analyse_elements):                      
                            # clean data
                            if(stats_colname == ''):  #data.shape[1] == 1):
                                data = data#.loc[:id_colname].to_csv(header=None, index=False).strip('\n').split('\n')
                                method = 'ora_df'
                                print(data)
                            else:
                                data = data[[id_colname, stats_colname]]#, data.columns[1]
                                method = 'ulm'
                            data[id_colname] = data[id_colname].apply(str.upper)
                            datasets.append({'datasetname': filename, 'data': data, 'method': method, 'datasetid': i})
    
    return datasets




def get_testdata(w_inputformat, datarootpath):
    """Read and process Testdata"""
    data = ''
    match w_inputformat:
        case UiVal.KINASES:
            data = ['Tex2_S265', 'Abl2_S822', 'Abl2_S819', 'Prkd3_S213', 'Tfeb_S108', 'Prkd3_S216', 'Foxk2_S389', 'Mtor_T2473', 'Peak1_S1203', 'Peak1_S1206', 'Lsr_S591', 'Eps15l1_S253', 'Lsr_S588', 'Dennd6a_S16', 'Unc5b_S531', 'Nfia_S310', 'Mapk6_S683', 'Mapk6_S687', 'Unc5b_S528', 'Tns3_Y773', 'Tfeb_S121', 'Mybbp1a_T1260', 'Washc2_S1048', 'Washc2_S1049', 'Tfeb_S113', 'Clip1_S194', 'Iws1_S667', 'Larp1_S743', 'Znf608_S626', 'Trmt10c_S79', 'Gpatch8_S733', 'Gpatch8_S735', 'Fnbp1_S299', 'Nfia_S303', 'Ankrd17_S1692', 'Rhbdf1_T180', 'Pnisr_S321', 'Trmt10c_S85', 'Micall1_S494', 'Micall1_S496', 'Ankrd17_S1696', 'Top2b_S1448', 'Rabep1_S410', 'Znf608_T635', 'Tnks1bp1_S889', 'Mrpl11_S160', 'Nufip2_T221', 'Morc2a_S737', 'Rgl3_S568', 'Rgl3_S572', 'Ercc5_S572', 'Lsm11_S18', 'Arhgap12_S211', 'Rbm25_S678', 'Samd4b_S585', 'Sipa1l1_S211', 'Supt5h_T661', 'Srrm2_S2691', 'Srrm2_T2689', 'Plec_S4396', 'Iws1_S666', 'Pnn_S698', 'Pnn_S700', 'Lmna_T19', 'Phldb2_S71']
            elementtype = 'kinases'
            elements = "'Tex2_S265', 'Abl2_S822', 'Abl2_S819', 'Prkd3_S213', 'Tfeb_S108', 'Prkd3_S216', 'Foxk2_S389', 'Mtor_T2473', 'Peak1_S1203', 'Peak1_S1206', 'Lsr_S591', 'Eps15l1_S253', 'Lsr_S588', 'Dennd6a_S16', 'Unc5b_S531', 'Nfia_S310', 'Mapk6_S683', 'Mapk6_S687', 'Unc5b_S528', 'Tns3_Y773', 'Tfeb_S121', 'Mybbp1a_T1260', 'Washc2_S1048', 'Washc2_S1049', 'Tfeb_S113', 'Clip1_S194', 'Iws1_S667', 'Larp1_S743', 'Znf608_S626', 'Trmt10c_S79', 'Gpatch8_S733', 'Gpatch8_S735', 'Fnbp1_S299', 'Nfia_S303', 'Ankrd17_S1692', 'Rhbdf1_T180', 'Pnisr_S321', 'Trmt10c_S85', 'Micall1_S494', 'Micall1_S496', 'Ankrd17_S1696', 'Top2b_S1448', 'Rabep1_S410', 'Znf608_T635', 'Tnks1bp1_S889', 'Mrpl11_S160', 'Nufip2_T221', 'Morc2a_S737', 'Rgl3_S568', 'Rgl3_S572', 'Ercc5_S572', 'Lsm11_S18', 'Arhgap12_S211', 'Rbm25_S678', 'Samd4b_S585', 'Sipa1l1_S211', 'Supt5h_T661', 'Srrm2_S2691', 'Srrm2_T2689', 'Plec_S4396', 'Iws1_S666', 'Pnn_S698', 'Pnn_S700', 'Lmna_T19', 'Phldb2_S71'"
            st.markdown(f"**The following {elementtype} are used:**")
            st.write(elements)
        case UiVal.GENES:
            data = ['KIAA0907', 'KDM5A', 'CDC25A', 'EGR1', 'GADD45B', 'RELB', 'TERF2IP', 'SMNDC1', 'TICAM1', 'NFKB2', 'RGS2', 'NCOA3', 'ICAM1', 'TEX10', 'CNOT4', 'ARID4B', 'CLPX', 'CHIC2', 'CXCL2', 'FBXO11', 'MTF2', 'CDK2', 'DNTTIP2', 'GADD45A', 'GOLT1B', 'POLR2K', 'NFKBIE', 'GABPB1', 'ECD', 'PHKG2', 'RAD9A', 'NET1', 'KIAA0753', 'EZH2', 'NRAS', 'ATP6V0B', 'CDK7', 'CCNH', 'SENP6', 'TIPARP', 'FOS', 'ARPP19', 'TFAP2A', 'KDM5B', 'NPC1', 'TP53BP2', 'NUSAP1', 'SCCPDH', 'KIF20A', 'FZD7', 'USP22', 'PIP4K2B', 'CRYZ', 'GNB5', 'EIF4EBP1', 'PHGDH', 'RRAGA', 'SLC25A46', 'RPA1', 'HADH', 'DAG1', 'RPIA', 'P4HA2', 'MACF1', 'TMEM97', 'MPZL1', 'PSMG1', 'PLK1', 'SLC37A4', 'GLRX', 'CBR3', 'PRSS23', 'NUDCD3', 'CDC20', 'KIAA0528', 'NIPSNAP1', 'TRAM2', 'STUB1', 'DERA', 'MTHFD2', 'BLVRA', 'IARS2', 'LIPA', 'PGM1', 'CNDP2', 'BNIP3', 'CTSL1', 'CDC25B', 'HSPA8', 'EPRS', 'PAX8', 'SACM1L', 'HOXA5', 'TLE1', 'PYGL', 'TUBB6', 'LOXL1']
            elementtype = 'genes'
            elements = "'KIAA0907', 'KDM5A', 'CDC25A', 'EGR1', 'GADD45B', 'RELB', 'TERF2IP', 'SMNDC1', 'TICAM1', 'NFKB2', 'RGS2', 'NCOA3', 'ICAM1', 'TEX10', 'CNOT4', 'ARID4B', 'CLPX', 'CHIC2', 'CXCL2', 'FBXO11', 'MTF2', 'CDK2', 'DNTTIP2', 'GADD45A', 'GOLT1B', 'POLR2K', 'NFKBIE', 'GABPB1', 'ECD', 'PHKG2', 'RAD9A', 'NET1', 'KIAA0753', 'EZH2', 'NRAS', 'ATP6V0B', 'CDK7', 'CCNH', 'SENP6', 'TIPARP', 'FOS', 'ARPP19', 'TFAP2A', 'KDM5B', 'NPC1', 'TP53BP2', 'NUSAP1', 'SCCPDH', 'KIF20A', 'FZD7', 'USP22', 'PIP4K2B', 'CRYZ', 'GNB5', 'EIF4EBP1', 'PHGDH', 'RRAGA', 'SLC25A46', 'RPA1', 'HADH', 'DAG1', 'RPIA', 'P4HA2', 'MACF1', 'TMEM97', 'MPZL1', 'PSMG1', 'PLK1', 'SLC37A4', 'GLRX', 'CBR3', 'PRSS23', 'NUDCD3', 'CDC20', 'KIAA0528', 'NIPSNAP1', 'TRAM2', 'STUB1', 'DERA', 'MTHFD2', 'BLVRA', 'IARS2', 'LIPA', 'PGM1', 'CNDP2', 'BNIP3', 'CTSL1', 'CDC25B', 'HSPA8', 'EPRS', 'PAX8', 'SACM1L', 'HOXA5', 'TLE1', 'PYGL', 'TUBB6', 'LOXL1'"
            st.markdown(f"**The following {elementtype} are used:**")
            st.write(elements)
        case UiVal.MATRIX:
            data = pd.read_csv(datarootpath+'/differential_stats.csv')
            datasets= {'datasetname': 'differential_stats', 'data': data, 'method': 'ulm', 'datasetid': 0}

        case UiVal.H5AD: 
            process_sc()
        case UiVal.EXCEL:
            data = pd.read_excel(datarootpath +'/dePhosRes.xlsx')
            data = data.loc[:,['site','t_statistic']]
            data['site'] = data['site'].apply(lambda x: x.upper())
            data[id_colname] = data[id_colname].apply(str.upper)
            datasets= {'datasetname': 'smartphos', 'data': data, 'method': 'ulm', 'datasetid': 0}
    return datasets


