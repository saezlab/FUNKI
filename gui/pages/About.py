#from copy import deepcopy
#from io import StringIO 
#from IPython.display import display, Markdown   #, Latex # to display Markdown in code chunk
#import matplotlib.pyplot as plt
#import matplotlib.image as mpimg
#from pdf2image import convert_from_path
#from PIL import Image
#from pdf2image import convert_from_path
#import subprocess
#subprocess.run()
#process = subprocess.run(['poetry run streamlit run', '/Users/hanna/Documents/projects/SGUI/CTLA4/v00/scripts/python/Main.py'], 
#                         stdout=subprocess.PIPE, 
#                         universal_newlines=True)
#process
# Read png
#im = Image.open("./logo_funki.png") # read the image, provide the correct path
#im.show() # display the image 
import streamlit as st
from os import makedirs, path
from copy import deepcopy
from IPython.display import display, Markdown   #, Latex # to display Markdown in code chunk                          # sys, numba, logging, random, re, dill, logging.config
from standard_workflows import * 
from pdf2image import convert_from_path
if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from utils import utilities as util
from streamlit_extras.switch_page_button import switch_page
from st_pages import Page, show_pages, hide_pages
from streamlit_elements import elements, mui

util.add_logo()

st.markdown('# About')
st.markdown("FUNKI provides methods for functional omics analysis. It uses prior knowledge from the database Omnipath and uses the statistical methods from the tool Decoupler to find footprints in the data, like transcription factors and pathways.\nFUNKI is still under development so that new features will be added step by step.")
dc = st.button("Decoupler")
if dc:
    switch_page('decoupler')
#drug = st.button("Drugst.One")
#if drug:
#    switch_page('Drugst.One')




# mui buttons but how can they be clicked??
# with elements("multiple_children"):
#     mui.Button(
#         mui.icon.DoubleArrow,
#         "Button"
#     )
#     mui.Button(
#         mui.icon.Arrow,
#         "Button"
#     )
# shoudl work but doesn't
#st.button('click',on_click=switch_page('Drugst.One'),args=('click'))