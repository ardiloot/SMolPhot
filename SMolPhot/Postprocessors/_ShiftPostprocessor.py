import logging
from _Common import PostprocessorBase
from SMolPhot.Components.BaseClasses import Param

__all__ = ["ShiftPostprocessor"]

class ShiftPostprocessor(PostprocessorBase):
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Postprocessor.ShiftPostprocessor")
                
        paramsThis = [Param("_dx", 0e-9, \
                        friendlyName = "x", \
                        guiSiPrefix = True, \
                        guiSuffix = "m", \
                        guiStep = 1e-9),
                      
                      Param("_dy", 0e-9, \
                        friendlyName = "y", \
                        guiSiPrefix = True, \
                        guiSuffix = "m", \
                        guiStep = 1e-9),
                      
                      Param("_dz", 0e-9, \
                        friendlyName = "z", \
                        guiSiPrefix = True, \
                        guiSuffix = "m", \
                        guiStep = 1e-9),
                      
                        Param("_dI", 0.0, \
                        friendlyName = "I", \
                        guiSuffix = " ph", \
                        guiStep = 10.0),
                      
                        Param("_multI", 1.0, \
                        friendlyName = "Multiply I", \
                        guiLimits = (1e-6, 100), 
                        guiStep = 0.01)]
        
        PostprocessorBase.__init__(self, paramsThis + params, **kwargs)
    
    def Apply(self, fseries, psf, groundTruth, locs):
        if not self._enabled:
            return locs

        for loc in locs:
            (x, y, z), photons = loc.coord, loc.photons
            newX = x + self._dx
            newY = y + self._dy
            newZ = z + self._dz if z is not None else None
            loc.coord = (newX, newY, newZ)
            
            if photons is not None:
                newPhotons = self._multI * (photons + self._dI)
                loc.photons = newPhotons
        return locs
        
    @property
    def name(self):
        return "Shift"
    