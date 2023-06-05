import webbrowser, glob              # os, json, yaml, base64, re, sys, numba, logging, random, dill, logging.config
from os import makedirs, path
import streamlit as st

from standard_workflows import * 
if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from utils import utilities as util


### BasePage ###
util.add_logo()


from string import Template
net = {"nodes": [{"id":"PTEN"}, {"id":"custom123"}], "edges": [{"from":"PTEN", "to":"custom123"}]}
html = """
<html>
<head>
   <script src="https://cdn.drugst.one/latest/drugstone.js"></script>
   <link rel="stylesheet" href="https://cdn.drugst.one/latest/styles.css">
   <style>
            h1 {
                font-family: "Arial";
            }
            .css-z5fcl4 {
                 padding: 2rem 2rem;
            }
        </style>
</head>
<body>
    <h1>Drugst.One Network Analysis</h1>
  <drugst-one
    config="{}"
    groups="{}"
    network="$net"> 
   </drugst-one>
</body>
</html>
"""
html = Template(html).safe_substitute(net=net)
st.components.v1.html(html, width=1000, height=1000, scrolling=False)





