import logging
import numpy as np
from scipy import interpolate
from _Base import PreprocessorBase
from SMolPhot.Components.BaseClasses import Param

__all__ = ["BackgroundSubractor"]

#===============================================================================
# Methods
#===============================================================================

def Bin2dData(data, bins, padding = True):
    data = np.asarray(data)
    bins = np.asarray(bins)
    inputShape = np.asarray(data.shape)
    elementsPerBin = np.ceil(inputShape.astype(float) / bins).astype(int)
    toPad = bins * elementsPerBin - inputShape
    
    if not (toPad == 0).all():
        if not padding:
            raise ValueError("Enable padding to bin this data")
        
        padBeginning = toPad / 2
        padEnd = toPad - padBeginning
        padWidth = ((padBeginning[0], padEnd[0]),(padBeginning[1], padEnd[1]))
        dataPadded = np.pad(data, padWidth, "constant", constant_values=(np.NaN, np.NaN))
    else:
        dataPadded = data
   
    res = dataPadded.reshape([bins[0], elementsPerBin[0], bins[1], elementsPerBin[1]])
    return res

#===============================================================================
# BackgroundSubractor
#===============================================================================

class BackgroundSubractor(PreprocessorBase):   

    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Preprocessors.BackgroundSubractor")
        self.logger.warn("Currently works only if it is the first preprocessor")
        paramsThis = [Param("_binsX", 4, friendlyName = "Bins X", guiStep = 1, guiLimits=(4, 1000)),
                      Param("_binsY", 4, friendlyName = "Bins Y", guiStep = 1, guiLimits=(4, 1000)),
                      Param("_spanT", 1, friendlyName = "Span T", guiStep = 1, guiLimits=(1, 10000)),
                      Param("_precentile", 50, friendlyName = "Percentile", guiStep = 1, guiLimits=(1, 99)),
                      Param("_smoothing", 0.0, friendlyName = "Smoothing")]
        #, friendlyName = "Tolerance XY", guiSiPrefix = True, guiSuffix = "m", guiStep = 10e-9
        PreprocessorBase.__init__(self, paramsThis + params, **kwargs)
    
    def Apply(self, frame):
        if not self.enabled:
            return
        self.logger.debug("Apply")
        
        bins = (self._binsX, self._binsY)

        # Find bin's central coordinates
        coordsMatX = np.nanmean(Bin2dData(frame.coords[0], bins), axis = (1, 3))
        coordsMatY = np.nanmean(Bin2dData(frame.coords[1], bins), axis = (1, 3))

        # Calc percentile over self._spanT frames
        # TODO: performance issue if _spanT large
        framesFrom = max(0, frame.nr - int(self._spanT / 2))
        framesTo = min(framesFrom + self._spanT, len(frame._fseries))
        
        binnedData = None
        for frameNr in range(framesFrom, framesTo):
            # TODO: if not first preprocessor
            dataF = frame._fseries.GetOriginalFrame(frameNr).data
            _binnedData = Bin2dData(dataF, bins)
            
            if binnedData is None:
                binnedData = _binnedData
            else:
                binnedData = np.append(binnedData, _binnedData, axis = 1)

        percentileMat = np.nanpercentile(binnedData, self._precentile, axis = (1, 3))
        
        # Interpolate plane
        #x, y, z = 1e9 * coordsMatX.ravel(), 1e9 * coordsMatY.ravel(), percentileMat.ravel()
        #tck = interpolate.bisplrep(x, y, z, s = self._smoothing)
        #interpolatedBg = interpolate.bisplev(1e9 * xc, 1e9 * yc, tck)
        
        xc, yc = frame.coords
        interp = interpolate.RectBivariateSpline(coordsMatX[:, 0], coordsMatY[0, :], percentileMat, s = self._smoothing)
        interpolatedBg = interp.ev(xc, yc)

        # Result
        frame.data -= interpolatedBg
        frame._interpolatedBg = interpolatedBg # TODO: 

        self.logger.debug("Apply done")
        return frame
    
    @property
    def name(self):
        return "Background subtractor"

if __name__ == "__main__":
    pass