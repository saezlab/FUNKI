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
#import utilities as util
import pandas as pd
import decoupler as dc
#import openpyxl
import plotly.express as px
import numpy as np
import kaleido, subprocess, glob
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

    priorKnwldg = sc_funcs.deep_get(ap['dataset_params'], (sc_funcs.getpath(ap['dataset_params'], 'priorKnowledge', search_value=False))).keys()
    organism = list(ap['dataset_params'].keys())[0]

    for resource_type in priorKnwldg:
        resources = sc_funcs.deep_get(ap['dataset_params'], sc_funcs.getpath(ap['dataset_params'], resource_type, search_value=False)).keys()
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
            targetfigurepaths=[]
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
                    
                    #st.write(result)
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
                        df = df.sort_values(by='t', axis=0, ascending=False, key = abs)
                        return df
                    df = decouple_result_totable(result)
                    
                    dc.plot_barplot(result[f'{method}_estimate'], result[f'{method}_estimate'].index[0], top=25, vertical=False, return_fig = False, save = figpath)
                    st.write(data)
                    data.index = data[data.columns[0]]
                    targets_to_inspect = df.index
                    if resource != 'ksn':
                        if len(targets_to_inspect) >= 20:
                            targets_to_inspect = targets_to_inspect[:20]
                        for source_name in targets_to_inspect:  #[:9]: #iloc[:,0][:5]:
                            targetfigurepath = f'{resource}_{source_name}.png'
                            dc.plot_targets(data, stat=data.columns[1], source_name=source_name, net=net, top=15, return_fig = False, save = targetfigurepath)
                            targetfigurepaths = targetfigurepaths + [targetfigurepath]
                    else:
                        targetfigurepaths = []
                except ValueError as ve:
                    try:
                        result = dc.decouple(d, net, methods=[method], min_n=2)
                        st.warning("Warning: All sources with at least two targets were taken into consideration. The default would be to have at least five targets per source. To go with the default you would need to add more genes.")
                    except ValueError as ve:
                        st.error("Error: There aren't any sources with the minimum of five targets. Please provide more genes.")
            util.add_results(df, figpath, title, f"{title}{dataset['datasetid']}", targetfigurepaths)
            


def get_data(w_inputformat)->list[dict] :
    datasets = list()
    def increment_analysiscount():
        st.session_state.analysiscount += 1
    match w_inputformat:
        case UiVal.H5AD: 
            st.caption("Please provide a .h5ad file that contains the log transformed data in its raw attribute. If you have single cell data that doesn't meet this condition, please contact Hanna.")
            uploaded_files = st.file_uploader("Choose a file", accept_multiple_files=True, type = ['.xlsx', '.h5ad']) 
            #st.write(uploaded_files)
            if len(uploaded_files) >= 1:
                analysispath = st.session_state.ap['proj_params']['paths']['analysis_path']
                st.session_state.ap['proj_params']['paths']['data_root_path'] = ''
                import scanpy
                data = scanpy.read_h5ad(uploaded_files[0])
                #st.write(data)
                proj_id = st.session_state.ap['proj_params']['proj_id']
                ap_dir =  os.path.join(analysispath, proj_id, 'v01', 'analysis')
                ap_path = os.path.join(ap_dir, 'analysis_params.yaml')
                import yaml
                if not os.path.exists(ap_dir):
                    os.makedirs(ap_dir)
                with open(ap_path, 'w+') as file: 
                    yaml.dump(st.session_state.ap, file, sort_keys=False)
 
                with open(ap_path) as stream:
                    apsfile = yaml.load(stream, Loader=yaml.BaseLoader)
                st.write(apsfile)

                dc_dataset = sc_classes.Analysis.new_dataset(sc_classes.Baseanalysis, dcu.Decoupler) 
                analysis = sc_classes.Analysis(datasets=[
                                                ('01', 'SingleCellRNAseq', 'human', dc_dataset)
                                                ], params_path = ap_dir)
                analysis.datasets[0].data = data
                st.write(analysis.datasets[0].data)
                scl.analysis = analysis
                scl.get_acts()
                scl.plot_umap()
                #scl.plot_violin('seurat_clusters')
                scl.get_mean_acts()
                scl.plot_mean_acts()
                st.write(analysis.datasets[0].get_paths())
                import zipfile
                from zipfile import ZipFile

                def zipDirectory(path, zippedFileName):
                    # traversing over all the files in a directory
                    for folderName, subfolders, filenames in os.walk(path):
                        for file in filenames:
                            # adding the current file to be zipped
                            zippedFileName.write(os.path.join(folderName, file))
                    print('Zipped the files!')
                # calling the constructor to create an object of the ZipFile class.

                filename_zip = f'{proj_id}.zip'
                zippedFileName = zipfile.ZipFile(filename_zip, 'w', zipfile.ZIP_DEFLATED)

                # calling the zipDirectory function by passing the path of the current directory and zipped file name.
                zipDirectory(f'./{proj_id}', zippedFileName)
                zippedFileName.close()

                if os.path.exists(filename_zip,):
                    with open(filename_zip, 'rb') as fp:
                        btn = st.download_button(
                            label='Download ZIP',
                            data=fp,
                            file_name=filename_zip,
                            mime='application/zip'
                        )       
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
                                st.write('You can find the results in the "Results 1st dataset" tab above.')
                    return datasets
                            
                datasets = read_and_display_files(uploaded_files)
    return datasets




def get_testdata(w_inputformat, w_omicstype, datarootpath):
    """Read and process Testdata"""
    data = ''
    datasets = list()
    match w_inputformat:
        case UiVal.KINASES:
            data = ["TEX2_S265", "ABL2_S822", "ABL2_S819", "PRKD3_S213", "TFEB_S108", "PRKD3_S216", "FOXK2_S389", "MTOR_T2473", "PEAK1_S1203", "PEAK1_S1206", "LSR_S591", "EPS15L1_S253", "LSR_S588", "DENND6A_S16", "UNC5B_S531", "NFIA_S310", "MAPK6_S683", "MAPK6_S687", "UNC5B_S528", "TNS3_Y773", "TFEB_S121", "MYBBP1A_T1260", "WASHC2_S1048", "WASHC2_S1049", "TFEB_S113", "CLIP1_S194", "IWS1_S667", "LARP1_S743", "ZNF608_S626", "TRMT10C_S79", "GPATCH8_S733", "GPATCH8_S735", "FNBP1_S299", "NFIA_S303", "ANKRD17_S1692", "RHBDF1_T180", "PNISR_S321", "TRMT10C_S85", "MICALL1_S494", "MICALL1_S496", "ANKRD17_S1696", "TOP2B_S1448", "RABEP1_S410", "ZNF608_T635", "TNKS1BP1_S889", "MRPL11_S160", "NUFIP2_T221", "MORC2A_S737", "RGL3_S568", "RGL3_S572", "ERCC5_S572", "LSM11_S18", "ARHGAP12_S211", "RBM25_S678", "SAMD4B_S585", "SIPA1L1_S211", "SUPT5H_T661", "SRRM2_S2691", "SRRM2_T2689", "PLEC_S4396", "IWS1_S666", "PNN_S698", "PNN_S700", "LMNA_T19", "PHLDB2_S71", "PHLDB2_S73", "LARP1_S751", "USP6NL_S544", "USP6NL_S547", "BCL9L_S118", "MTOR_S2478", "MTOR_S2481", "FARP1_S376", "SMARCA4_S695", "FARP1_S373", "ACACA_S25", "ARHGEF17_T695", "DCP2_S246", "DCP2_S247", "DCP2_S249", "AFF4_S485", "AFF4_S486", "RGL3_S577", "OTUD4_S1005", "ARID4A_T1148", "SMARCC1_S25", "SMARCC1_Y30", "ARHGEF5_S593", "AFF4_S482", "COBLL1_S375", "FZR1_S146", "SOS1_S1251", "SOS1_T1249", "SNX2_T104", "SLAIN2_S59", "ECT2_S376", "CROCC_S1483", "TNS2_S825", "ERCC5_S571", "WRNIP1_S91", "WRNIP1_S92", "LRRC8A_S217", "LRRC8A_T215", "NFIC_S343", "MARCKS_S138", "RAI14_T297", "DPYSL2_T514", "UGDH_T474", "OSBPL8_S807", "OSBPL8_S808", "NFIC_S284", "SEC61B_S17", "ARID4B_S675", "CTTNBP2NL_S562", "ARHGEF17_T692", "CEP131_T473", "ANKS1_S665", "CELSR1_S2776", "AKAP8_S320", "AKAP8_S325", "KAT7_T87", "RRAGC_S2", "RABEP1_S407", "VIM_S2", "COBLL1_S374", "ARHGAP12_Y241", "DNMT1_S150", "HIRIP3_S396", "ETL4_S359", "ACACA_S23", "EEPD1_S31", "SLC12A6_S1032", "TRMT10C_S76", "BCL9L_S116", "SMARCC1_T815", "CEP350_S1256", "JUNB_S256", "GAS2L1_S482", "ATN1_S101", "CXADR_T334", "TNS2_S820", "PWP2_T891", "JMY_S130", "JMY_S134", "RERE_S594", "GAS2L1_S489", "AEBP2_S203", "CENPU_S182", "OSBPL8_S810", "DPYSL2_S522", "SAFB_S626", "FNBP1L_S501", "JUNB_T252", "ARID4A_S1149", "SLC7A6OS_S299", "BAIAP2L1_S413", "RFC1_S244", "PAK2_S152", "BOP1_S14", "SRSF10_S129", "WWC2_S665", "ATXN2_S702", "KIAA1522_S883", "KIAA1522_S886", "CGNL1_S284", "AGPS_S52", "MARCKS_S163", "MKI67_S498", "MKI67_S503", "PLEKHG3_S736", "PLEKHG3_S737", "GPATCH8_S1039", "ECT2_T373", "F11R_S285", "SPAG5_S12", "SPAG5_S14", "BOP1_S11", "ATRX_S1318", "CTTNBP2NL_S567", "MPRIP_T623", "SCRIB_S1284", "PRPF4B_S33", "CENPU_S186", "COBLL1_S396", "RCC2_S48", "RCC2_S49", "XRN2_S487", "SHROOM1_S49", "SHROOM1_T47", "CLCC1_S429", "CLCC1_S433", "CELSR1_S2779", "AGPS_S57", "OSTF1_S203", "TMEM50A_S8", "RCOR3_S227", "PHLDB2_S329", "ZYX_S336", "MORC2A_S741", "PRRC2A_S342", "PRRC2A_S350", "PNISR_S311", "PNISR_S313", "PHF23_S74", "PHF23_S79", "MYBBP1A_T1256", "DNMT1_S140", "MAP2_S1783", "ARFGEF2_S218", "COBLL1_S384", "MICALL2_S655", "TJP1_S131", "ITGB4_S1451", "ITGB4_S1454", "KRT8_S43", "GTPBP1_S25", "CWC22_S60", "TRIM28_S473", "CCM2_S417", "LRRC47_S293", "COBLL1_S366", "CMTM4_T208", "GTPBP1_S8", "SRRM2_S980", "SGTA_Y163", "FNBP1_S302", "ATN1_S103", "MARCKS_T143", "CUL4B_S110", "SAMD4B_S593", "FAM83H_S937", "TNKS1BP1_S1063", "GCFC2_S16", "ATXN2_S697", "MCM2_S5", "CMTR1_S52", "MYO18A_S157", "LARP1_S299", "RALGAPA2_S375", "RCOR3_S212", "DBNDD2_T137", "PCM1_S1777", "NFIC_S305", "LARP1_S302", "MACROH2A1_S170", "WWC2_S662", "RHBDF1_T183", "HELLS_S495", "CMTR1_S50", "NOLC1_T594", "PDCD4_S457", "OSBPL11_S179", "PA2G4_T11", "ATRX_T740", "LSM11_S15", "CROCC_S1479", "LMNA_S633", "KANSL3_S536", "KANSL3_S540", "ETV3_S250", "ETV3_S245", "PHACTR4_S148", "SEPTIN9_S85", "UHRF1_S91", "ARFGEF2_S227", "VCPIP1_S993", "FNBP1_S296", "CHD2_S156", "CHD2_S165", "CHD2_T164", "DPYSL2_S518", "PTMS_S5", "GIGYF2_S161", "GZF1_S193", "GZF1_T196", "CHD2_S168", "TOX4_S182", "PDS5B_S1176", "CXADR_S332"]
            elementtype = 'kinases'
            elements = "TEX2_S265, ABL2_S822, ABL2_S819, PRKD3_S213, TFEB_S108, PRKD3_S216, FOXK2_S389, MTOR_T2473, PEAK1_S1203, PEAK1_S1206, LSR_S591, EPS15L1_S253, LSR_S588, DENND6A_S16, UNC5B_S531, NFIA_S310, MAPK6_S683, MAPK6_S687, UNC5B_S528, TNS3_Y773, TFEB_S121, MYBBP1A_T1260, WASHC2_S1048, WASHC2_S1049, TFEB_S113, CLIP1_S194, IWS1_S667, LARP1_S743, ZNF608_S626, TRMT10C_S79, GPATCH8_S733, GPATCH8_S735, FNBP1_S299, NFIA_S303, ANKRD17_S1692, RHBDF1_T180, PNISR_S321, TRMT10C_S85, MICALL1_S494, MICALL1_S496, ANKRD17_S1696, TOP2B_S1448, RABEP1_S410, ZNF608_T635, TNKS1BP1_S889, MRPL11_S160, NUFIP2_T221, MORC2A_S737, RGL3_S568, RGL3_S572, ERCC5_S572, LSM11_S18, ARHGAP12_S211, RBM25_S678, SAMD4B_S585, SIPA1L1_S211, SUPT5H_T661, SRRM2_S2691, SRRM2_T2689, PLEC_S4396, IWS1_S666, PNN_S698, PNN_S700, LMNA_T19, PHLDB2_S71, PHLDB2_S73, LARP1_S751, USP6NL_S544, USP6NL_S547, BCL9L_S118, MTOR_S2478, MTOR_S2481, FARP1_S376, SMARCA4_S695, FARP1_S373, ACACA_S25, ARHGEF17_T695, DCP2_S246, DCP2_S247, DCP2_S249, AFF4_S485, AFF4_S486, RGL3_S577, OTUD4_S1005, ARID4A_T1148, SMARCC1_S25, SMARCC1_Y30, ARHGEF5_S593, AFF4_S482, COBLL1_S375, FZR1_S146, SOS1_S1251, SOS1_T1249, SNX2_T104, SLAIN2_S59, ECT2_S376, CROCC_S1483, TNS2_S825, ERCC5_S571, WRNIP1_S91, WRNIP1_S92, LRRC8A_S217, LRRC8A_T215, NFIC_S343, MARCKS_S138, RAI14_T297, DPYSL2_T514, UGDH_T474, OSBPL8_S807, OSBPL8_S808, NFIC_S284, SEC61B_S17, ARID4B_S675, CTTNBP2NL_S562, ARHGEF17_T692, CEP131_T473, ANKS1_S665, CELSR1_S2776, AKAP8_S320, AKAP8_S325, KAT7_T87, RRAGC_S2, RABEP1_S407, VIM_S2, COBLL1_S374, ARHGAP12_Y241, DNMT1_S150, HIRIP3_S396, ETL4_S359, ACACA_S23, EEPD1_S31, SLC12A6_S1032, TRMT10C_S76, BCL9L_S116, SMARCC1_T815, CEP350_S1256, JUNB_S256, GAS2L1_S482, ATN1_S101, CXADR_T334, TNS2_S820, PWP2_T891, JMY_S130, JMY_S134, RERE_S594, GAS2L1_S489, AEBP2_S203, CENPU_S182, OSBPL8_S810, DPYSL2_S522, SAFB_S626, FNBP1L_S501, JUNB_T252, ARID4A_S1149, SLC7A6OS_S299, BAIAP2L1_S413, RFC1_S244, PAK2_S152, BOP1_S14, SRSF10_S129, WWC2_S665, ATXN2_S702, KIAA1522_S883, KIAA1522_S886, CGNL1_S284, AGPS_S52, MARCKS_S163, MKI67_S498, MKI67_S503, PLEKHG3_S736, PLEKHG3_S737, GPATCH8_S1039, ECT2_T373, F11R_S285, SPAG5_S12, SPAG5_S14, BOP1_S11, ATRX_S1318, CTTNBP2NL_S567, MPRIP_T623, SCRIB_S1284, PRPF4B_S33, CENPU_S186, COBLL1_S396, RCC2_S48, RCC2_S49, XRN2_S487, SHROOM1_S49, SHROOM1_T47, CLCC1_S429, CLCC1_S433, CELSR1_S2779, AGPS_S57, OSTF1_S203, TMEM50A_S8, RCOR3_S227, PHLDB2_S329, ZYX_S336, MORC2A_S741, PRRC2A_S342, PRRC2A_S350, PNISR_S311, PNISR_S313, PHF23_S74, PHF23_S79, MYBBP1A_T1256, DNMT1_S140, MAP2_S1783, ARFGEF2_S218, COBLL1_S384, MICALL2_S655, TJP1_S131, ITGB4_S1451, ITGB4_S1454, KRT8_S43, GTPBP1_S25, CWC22_S60, TRIM28_S473, CCM2_S417, LRRC47_S293, COBLL1_S366, CMTM4_T208, GTPBP1_S8, SRRM2_S980, SGTA_Y163, FNBP1_S302, ATN1_S103, MARCKS_T143, CUL4B_S110, SAMD4B_S593, FAM83H_S937, TNKS1BP1_S1063, GCFC2_S16, ATXN2_S697, MCM2_S5, CMTR1_S52, MYO18A_S157, LARP1_S299, RALGAPA2_S375, RCOR3_S212, DBNDD2_T137, PCM1_S1777, NFIC_S305, LARP1_S302, MACROH2A1_S170, WWC2_S662, RHBDF1_T183, HELLS_S495, CMTR1_S50, NOLC1_T594, PDCD4_S457, OSBPL11_S179, PA2G4_T11, ATRX_T740, LSM11_S15, CROCC_S1479, LMNA_S633, KANSL3_S536, KANSL3_S540, ETV3_S250, ETV3_S245, PHACTR4_S148, SEPTIN9_S85, UHRF1_S91, ARFGEF2_S227, VCPIP1_S993, FNBP1_S296, CHD2_S156, CHD2_S165, CHD2_T164, DPYSL2_S518, PTMS_S5, GIGYF2_S161, GZF1_S193, GZF1_T196, CHD2_S168, TOX4_S182, PDS5B_S1176, CXADR_S332"
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


