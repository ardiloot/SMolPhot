"""This package contains different preprocessors. Preprocessors are separated
into different files to ease parallel work. The aim of the preprocessor is
to prepare frame for a localizer (reduce noise, transformations, etc.)

"""

from _BackgroundSubtractor import *
from _GaussianFilter import *
#from _MedianFilter import *
#from _WaveletFilter import *
#from _DenoiseFilter import *