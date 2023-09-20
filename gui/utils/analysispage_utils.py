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




def walking_footprints():
    """ Adds a html component with an animation of walking mouse footprints"""

    import streamlit.components.v1 as components
    components.html("""
            <style>
                            
                #leftfoot, #rightfoot, #frontleftfoot, #frontrightfoot {
                position:absolute;
                display:block;
                border: solid 1px transparent;
                -webkit-backface-visibility:hidden;
                }
                .foot {
                position:absolute;
                width:29px;
                height:29px;
                }

                #rightfoot{
                left:146px;
                top: 184px;
                visibility: hidden;
                width: 20px;
                height: 42px;
                }

                #leftfoot{
                left:90px;
                top: 184px;
                visibility: hidden;
                width: 20px;
                height: 42px;
                }

                #frontrightfoot{
                left:136px;
                top: 164px;
                visibility: hidden;
                width: 20px;
                height: 20px;
                }

                #frontleftfoot{
                left:95px;
                top: 164px;
                visibility: hidden;
                width: 20px;
                height: 20px;
                }


        </style>
    <img src="https://i.ibb.co/8xtMhp2/straight-Left-Foot.jpg" id="leftfoot">
    <img src="https://i.ibb.co/1Xh6M15/rightbackfootstraight.jpg" id="rightfoot">
    <img src="https://i.ibb.co/mH90bGM/leftfrontfootstraight.jpg" id="frontleftfoot">
    <img src="https://i.ibb.co/GWQ0q0W/frontrightfootstraight.jpg" id="frontrightfoot">

    <script src="https://cdnjs.cloudflare.com/ajax/libs/gsap/latest/TweenMax.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/gsap/2.0.2/TimelineMax.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/2.1.3/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/gsap/latest/utils/Draggable.min.js"></script>
    
    <script>
        var $rightfoot = $("#rightfoot"),
            $leftfoot = $("#leftfoot"),
            $frontrightfoot = $("#frontrightfoot"),
            $frontleftfoot = $("#frontleftfoot");

        var tl = new TimelineMax({repeat:0}) /*-1 to repeat indefinitely*/

        var ease = SteppedEase.config(5);

        tl

        .to($rightfoot, 0.25, {autoAlpha:1,},1)
        .to($leftfoot, 0.25, {autoAlpha:1,},1)

        .to($frontrightfoot, 0.25, {autoAlpha:1,},1)
        .to($frontleftfoot, 0.25, {autoAlpha:1,},1)

        .to($leftfoot, 3, 
        {bezier:{ curviness: 1, values:[{x:0, y:0},{x:52, y:-25}, {x:97, y:-39}, {x:136, y:-54}, {x:172, y:-83}, {x:197, y:-117},{x:200, y:-163,}],
        autoRotate:90}, ease:ease},1.5)

        .to($rightfoot, 3, 
        {bezier:{ curviness: 1, values:[{x:0, y:0}, {x:44, y:-13}, {x:97, y:-28,}, {x:134, y:-49,}, {x:166, y:-83,}, {x:182, y:-126,}, {x:176, y:-163,}],
        autoRotate:100}, ease:ease},1.75)

        .to($frontleftfoot, 3, 
        {bezier:{ curviness: 1, values:[{x:0, y:0},{x:52, y:-25}, {x:97, y:-39}, {x:136, y:-54}, {x:172, y:-83}, {x:197, y:-117},{x:200, y:-163,}],
        autoRotate:90}, ease:ease},1)

        .to($frontrightfoot, 3, 
        {bezier:{ curviness: 1, values:[{x:0, y:0}, {x:44, y:-13}, {x:97, y:-28,}, {x:134, y:-49,}, {x:166, y:-83,}, {x:182, y:-126,}, {x:176, y:-163,}],
        autoRotate:100}, ease:ease},1.25)
    </script>
    <div>
        <img src="https://i.ibb.co/8xtMhp2/straight-Left-Foot.jpg" id="leftfoot">
        <img src="https://i.ibb.co/1Xh6M15/rightbackfootstraight.jpg" id="rightfoot">
        <img src="https://i.ibb.co/mH90bGM/leftfrontfootstraight.jpg" id="frontleftfoot">
        <img src="https://i.ibb.co/GWQ0q0W/frontrightfootstraight.jpg" id="frontrightfoot">
    </div>
    """)



def init_page():
    st.markdown("# Transcription Factor & Pathway Analysis")
    coldesc, colfootprints = st.columns(2)
    with coldesc:
        # heading + description
        st.markdown("**Analysis of bulk RNA sequencing data with the tool [Decoupler](https://decoupler-py.readthedocs.io/en/latest/index.html).**\n\n   \
Start your analysis by uploading and describing your data in the 'Analysis' tab.\n   \
- The transcription factor (TF) analysis is done based on the prior knowledge resource `CollecTri`.\n  \
- The pathway analysis is based on the resource `Progeny`.")
                    
    with colfootprints:
        walking_footprints()
    st.markdown("---")