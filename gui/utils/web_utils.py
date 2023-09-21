import streamlit as st
import os, base64, re, sys
from os import makedirs, path
from copy import deepcopy
from IPython.display import display, Markdown   #, Latex # to display Markdown in code chunk
import json, yaml                               # numba, logging, random, dill, logging.config
from standard_workflows import *  
sys.path.append('../')
from streamlit_extras.switch_page_button import switch_page
from st_pages import Page, show_pages, _hide_pages

""" Add HTML, CSS, JS """

def add_logo():
    """ logo and css
    """
    #st.image(im, width = 200)
    st.markdown(
        """
        <style>
            h1 {
                font-family: "Arial";
                padding: 0rem 0rem 1rem;
                margin-top: 0.8rem;
                font-size: 2.2em;
                font-weight: 530;
            }

            .css-1a1tcp.e1ewe7hr3{
                color: black;
            }

            .css-18ni7ap {
                background: #ebebeb;
            }
            
            .css-14xtw13.e13qjvis0::before{
                content: "FUNKI ";
                float:left;
            }

            /*[role=tablist]{
                background-color: yellow: #FAC710 = rgba(250, 199, 16, .7); blue: #0280A3 = rgba(2, 128, 163, .8)*/
                /*background-color: rgba(255, 215, 0, .7);      /* last number is transparancy*/      
            }*/

            [role=tab]{
                background-color: rgba(0, 0, 0, 0);      /* make fully transparent so that the background of the tabs doesn't overlap/interfere with the background of the tablist*/      
            }

            [role=tab] div p{
                font-weight: 550;
            }

            .css-z5fcl4 {
                 padding: 2rem 2rem;
            }

            [data-testid="column"] {
                /*margin: 35px;
                padding: 10px;*/
            }

            [data-testid="stSidebar"] {
                background-color: #3a3a3a;
            }

            [data-testid="stSidebar"] span {
                color: #FFFFFF;
            }

            [data-testid="stSidebarNav"] {
                background-image: url('https://picr.eu/images/2023/06/11/sWpkI.jpg'); 
                background-size: contain;
                background-repeat: no-repeat;
                padding-top: 120px;
                background-position: 0px 45px;
                color: white;
            }
            [data-testid="stSidebarNav"]::before {
                /*content: "Methods";*/
                color: white;
                margin-left: 20px;
                margin-top: 20px;
                font-size: 30px;
                position: relative;
                top: 100px;
            }

            .css-l3i8zm {
                font-weight: 900;
                color: #08046e; 
            }

            [data-testid="stMarkdownContainer"] p {
                font-weight: 400;
            }

            [data-baseweb="select"] div {
                background-color: lightblue;
                width = min-content;
            }

            footer {
                visibility: hidden;
            }

            footer:after {
                content:'   Made with Streamlit by Hanna Schumacher, Copyright 2023, Heidelberg University Hospital';
                visibility: visible;
                display: block;
                height: 50px;
                clear: both;
                color: white;
                background-color: lightslategray;
            }

            /* markdown --- line */
            hr{
               margin: 0px;     
            }
            
            /* gaps between elements on the page */
            .css-pmz2b6{
                gap: 1rem;  
            }

            .css-fblp2m{
                fill: white;
            }
   
        </style>
        """,
        unsafe_allow_html=True,
    )



def init_page():
    """ logo, sidebar
    """
    add_logo()
    show_pages([
        Page("Analysis.py", "Analysis", ":bar_chart:"), 
        Page("gui/pages/About.py","About", "🏠"),
        #Page("gui/pages/Drugst.One.py","Drugst.One", ":spider_web:"),
        Page("gui/pages/Citations.py", "Citations", ":book:")
    ])
    #hide_pages(['utilities'])
    with st.sidebar:
        st.caption('FUNKI is under active development.</br>For the latest features see the <a href="https://saezlabfunkidev.streamlit.app/">development version</a>. </br>The code is on <a href="https://github.com/saezlab/FUNKI">Github</a>.', unsafe_allow_html=True)


