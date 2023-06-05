import streamlit as st
if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from utils import utilities as util


util.init_page()
st.markdown('# Citations')

st.markdown('## Funki')

st.write('<link to ShinyFunki paper>')

st.markdown('## Decoupler')
st.write('Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov D., Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores R.O. and Saez-Rodriguez J. 2022. decoupleR: ensemble of computational methods to infer biological activities from omics data. Bioinformatics Advances. https://doi.org/10.1093/bioadv/vbac016')

st.markdown('**License**')

st.markdown('Footprint methods inside decoupler can be used for academic or commercial purposes.') #, except ``viper`` which holds a non-commercial license.') 

st.markdown('## Omnipath')
st.markdown('Collectri and Progeny are ressources derived from Omnipath via the tool Decoupler.') 
st.markdown('The data redistributed by OmniPath does not have a license, each original resource carries their own. `Here <https://omnipathdb.org/info>`_ one can find the license information of all the resources in OmniPath.')



