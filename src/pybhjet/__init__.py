
# enables to directly use eg.: from pybhjet import preprocess_component_output
# instead of from pybhjet.bhjet_plotting import preprocess_component_output

from .pybhjet import *  
from .bhjet_plotting import * 

# this leads to 3ml being imported with every pybhjet import
# from .pybhjet_3ml import *