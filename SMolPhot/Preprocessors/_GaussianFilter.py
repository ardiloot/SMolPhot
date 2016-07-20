import logging
from _Base import PreprocessorBase
from SMolPhot.Components.BaseClasses import Param
from scipy.ndimage import filters
import numpy as np

__all__ = ["GaussianFilterPreprocessor"]

class GaussianFilterPreprocessor(PreprocessorBase):
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Preprocessors.GaussianFilterPreprocessor")
        paramsThis = [Param("_sigma", 50e-9, friendlyName = "Sigma", \
                            guiSiPrefix = True, guiSuffix = "m", \
                            guiStep = 1e-9, guiLimits = (1e-9, 1.0))]
        PreprocessorBase.__init__(self, paramsThis + params, **kwargs)
    
    def Apply(self, frame):
        if not self.enabled:
            return
        self.logger.debug("Apply")
        sigma = self._sigma / frame._fseries.pixelSize
        frame.data = filters.gaussian_filter(frame.data, sigma)

        self.logger.debug("Apply done")
        
    @property
    def name(self):
        return "Gaussian filter"
    