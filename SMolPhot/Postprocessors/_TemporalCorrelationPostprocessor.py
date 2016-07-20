import logging
from _Common import PostprocessorBase
from SMolPhot.Components.BaseClasses import Param
from collections import OrderedDict

__all__ = ["TemporalCorrelationPostprocessor"]


class TemporalCorrelationPostprocessor(PostprocessorBase):
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Postprocessor.TemporalCorrelationPostprocessor")

        self._weightComboValues = OrderedDict([("No weights", lambda loc: 1.0),
                                               ("PSF goodness", lambda loc: loc.goodnessValue),
                                               ("Photons", lambda loc: loc.photons),
                                               ("stdAmp", lambda loc: 1.0 / loc.fitStdAmp),
                                               ("stdOffset", lambda loc: 1.0 / loc.fitStdOffset)])

        paramsThis = [Param("_xyTol",
                            80e-9,
                            friendlyName = "xy Tol",
                            guiSiPrefix = True,
                            guiSuffix = "m",
                            guiStep = 1e-9),

                      Param("_darkFrames",
                            0,
                            friendlyName="Dark frames",
                            guiStep=1),

                      Param("_weights",
                            "No weights",
                            friendlyName = "Weights",
                            guiType = "list",
                            guiValues=self._weightComboValues.keys()
                            )]

        PostprocessorBase.__init__(self, paramsThis + params, **kwargs)
    
    def Apply(self, fseries, psf, groundTruth, locs):
        if not self._enabled:
            return locs
        
        weightFun = self._weightComboValues[self._weights]
        locBuffer = []
        for loc in locs:
            if len(locBuffer) == 0:  # buffer is empty
                locBuffer.append(BufferBin(loc, weightFun))
            else:
                # Remove old bins
                locBufferNew = []
                for locBin in locBuffer:
                    if loc.frameNr - locBin.lastFrameNr > self._darkFrames + 1:
                        # large temporal gap: clear bin
                        locBin.SaveWeightedLocs()
                        continue
                    locBufferNew.append(locBin)
                locBuffer = locBufferNew
                
                # Check for matching bin
                locationAdded = False
                for locBin in locBuffer:
                    if (loc.x - locBin.x) ** 2 + \
                        (loc.y - locBin.y) ** 2 < self._xyTol ** 2:
                        # location fits into an existing bin
                        locBin.Add(loc)
                        locationAdded = True
                        break
                    
                # If no matching bins: create a new bin
                if not locationAdded:
                    locBuffer.append(BufferBin(loc, weightFun))  

        # Save weighted locs of the remaining bins
        for locBin in locBuffer:
            locBin.SaveWeightedLocs()

        return locs
        
    @property
    def name(self):
        return "Temporal correlation"


class BufferBin:

    def __init__(self, loc, weightFunc):
        self.locs = []
        self.lastFrameNr = None
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.weightFunc = weightFunc
        self.weightSum = 0.0
        
        self.Add(loc)

    def Add(self, loc):
        self.locs.append(loc)
        self.lastFrameNr = loc.frameNr
        weight = float(self.weightFunc(loc))
        
        # weighted average
        newWeightSum = self.weightSum + weight  # add new weight
        self.x = (self.weightSum * self.x + weight * loc.x) / newWeightSum
        self.y = (self.weightSum * self.y + weight * loc.y) / newWeightSum
        
        # Check if we have z-coordinate
        if loc.z is not None:
            self.z = (self.weightSum * self.z + weight * loc.z) / newWeightSum
        else:
            self.z = None
            
        self.weightSum = newWeightSum

    def SaveWeightedLocs(self):
        if len(self.locs) > 1:
            for loc in self.locs:
                loc.coord = self.x, self.y, self.z
