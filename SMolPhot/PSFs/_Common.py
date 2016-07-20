import copy
import numpy as np

from collections import defaultdict
from scipy import interpolate, optimize
from SMolPhot.Components.BaseClasses import ParamsBaseClass

__all__ = []

#===============================================================================
# MoleculeLocation
#===============================================================================

class MoleculeLoc():
        
    def __init__(self, coord, coordPixels = None, initial = None, frameNr = None, \
                 bestFitParams = None, **kwargs):
        
        self.coord = coord
        self.coordPixels = coordPixels
        self.initial = (None, None, None) if initial is None else initial
        self.frameNr = frameNr
        self.bestFitParams = bestFitParams
        self.bboxHWx, self.bboxHWy = 5, 5
        
        self.residual = kwargs.pop("residual", None)
        self.fitStds = kwargs.pop("fitStds", None)
        self.photons = kwargs.pop("photons", 0.0)
        self.amp = kwargs.pop("amp", None)
        self.sigmaX = kwargs.pop("sigmaX", None)
        self.sigmaY = kwargs.pop("sigmaY", None)
        self.offset = kwargs.pop("offset", None)
        
        self.fitStdAmp = kwargs.pop("fitStdAmp", None)
        self.fitStdX0 = kwargs.pop("fitStdX0", None)
        self.fitStdY0 = kwargs.pop("fitStdY0", None)
        self.fitStdZ0 = kwargs.pop("fitStdZ0", None)
        self.fitStdOffset = kwargs.pop("fitStdOffset", None)
        
        self.iteration = kwargs.pop("iteration", None)
        self.fitPixels = kwargs.pop("fitPixels", None)
        
        
        self.startPositionNr = None
        self.startPositionOffset = (None, None)  # XY distance from ground-truth
        self.distXY = np.inf  # XY distance from ground truth
        self.distZ = np.inf  # Z distance from ground truth
        self.minFitDistXY = np.inf  # minimal XY distance from nearest molecule in frame
        self.goodnessValue = None  # Precalculated goodness
        self.errorValue = None  # Precalculated error probability
        
    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]
    
    @property
    def z(self):
        return self.coord[2]
    
    @property
    def pixelX(self):
        return self.coordPixels[0]

    @property
    def pixelY(self):
        return self.coordPixels[1]
    
    @property
    def absz(self):
        return abs(self.coord[2])
    
    @property
    def ix(self):
        return self.initial[0]

    @property
    def iy(self):
        return self.initial[1]
    
    @property
    def iz(self):
        return self.initial[2]
        
    @property
    def lsError(self):
        return np.sum(self.residual ** 2.0)
    
    @property
    def fitLocStd(self):
        return np.sqrt(self.fitStdX0 ** 2.0 + self.fitStdY0 ** 2.0)
    
    @property
    def nrOfPixels(self):
        return self.residual.size
    
    @property
    def absoffset(self):
        return abs(self.offset)

#===============================================================================
# AxialCalibParam
#===============================================================================

class AxialCalibParam(object):
    
    def __init__(self, calibPoints, valueFunc, smoothing, useStds = False):
        self._smoothing = None
        self.useStds = useStds
        
        # Build temproary datastructure by fitPixels
        valuesTmp = defaultdict(lambda : defaultdict(list))
        for frame, locs in calibPoints:
            for loc in locs:
                value = valueFunc(loc)
                valuesTmp[loc.fitPixels][frame.z].append(value)
                
        # Find zs, means, stds, ...
        self._dataDict = {}
        for fitPixels, valuesByZ in valuesTmp.iteritems():
            zs = np.array(sorted(valuesByZ.keys()))
            values = map(np.array, [valuesByZ[z] for z in zs])
            means = np.array([np.mean(v) for v in values])
            stds = np.array([np.std(v) for v in values]) if useStds else 4e-8 * np.ones_like(zs)
            stdMax = stds.max()
            stds = np.where(stds == 0.0, 4e-8 if stdMax == 0.0 else stdMax, stds)
            
            self._dataDict[fitPixels] = {"zs": zs, "values": values, "means": means, "stds": stds}
        
        self._maxFitPixels = max(self.fitPixelsList)
        self._minFitPixels = min(self.fitPixelsList)
        
        self.SetSmoothing(smoothing)
        
    def SetSmoothing(self, smoothing):
        if smoothing == self._smoothing:
            return
        
        # Update interpolation
        self._smoothing = smoothing
        for _, data in self._dataDict.iteritems():
            intCoefs = interpolate.splrep(data["zs"], data["means"], \
                                          w = 1.0 / data["stds"], s = smoothing)
            data["interpCoefs"] = intCoefs
            
    def GetInterpolatedValue(self, z, fitPixels):
        fitPixelsActual = min(self._maxFitPixels, max(self._minFitPixels, fitPixels))
        res = splevFast(z, self._dataDict[fitPixelsActual]["interpCoefs"])
        return res
    
    def GetData(self, fitPixels, error = True):
        if not error:
            fitPixelsActual = min(self._maxFitPixels, max(self._minFitPixels, fitPixels))
        else:
            fitPixelsActual = fitPixels
        res = self._dataDict[fitPixelsActual]
        return res
    
    @property
    def fitPixelsList(self):
        return sorted(self._dataDict.keys())

#===============================================================================
# PsfBase
#===============================================================================

class PsfBase(ParamsBaseClass):
    """Baseclass for all PSF calculations.
    
    """
    
    def __init__(self, *args, **kwargs):
        self._hasCalibrationData = False
        self.SetAxialCalibPoints([])
        ParamsBaseClass.__init__(self, *args, **kwargs)

    @property
    def hasCalibrationData(self):
        return self._hasCalibrationData
  
    @property
    def name(self):
        raise NotImplementedError()
    
#===============================================================================
# Methods
#===============================================================================

def GetLocsFrameNrs(locs):
    return np.array([loc.frameNr for loc in locs], dtype = int)

def GetLocsCoordsArray(locs):
    res = np.zeros((len(locs), 3))
    for i, loc in enumerate(locs):
        res[i, 0] = loc.x
        res[i, 1] = loc.y
        res[i, 2] = loc.z if loc.z is not None else 0.0
    return res

def Gaussian2DSym((xs, ys), amp, x0, y0, sigma, offset):
    x, y = xs - x0, ys - y0
    res = amp * np.exp(-((x) ** 2.0) / (2.0 * sigma ** 2.0) - \
                            ((y) ** 2.0) / (2.0 * sigma ** 2.0)) + offset
    return res

def LocsToArrays(locs):
    nrs = np.zeros((len(locs),), dtype = int)
    coordsX = np.zeros((len(locs),), dtype = float)
    coordsY = np.zeros((len(locs),), dtype = float)
    coordsZ = np.zeros((len(locs),), dtype = float)
    photons = np.zeros((len(locs),), dtype = float)
    
    for i, loc in enumerate(locs):
        nrs[i] = loc.frameNr
        coordsX[i] = loc.x
        coordsY[i] = loc.y
        coordsZ[i] = 0.0 if loc.z is None else loc.z
        photons[i] = loc.photons
    
    return nrs, (coordsX, coordsY, coordsZ), photons

def DuplicateLocsList(locsToCopy):
    # Does shallow copy and replaces coord (may be modified turing postprocessing)
    res = []
    for original in locsToCopy:
        duplicate = copy.copy(original)
        originalZ = float(original.z) if original.z is not None else None
        duplicate.coord = (float(original.x), float(original.y), originalZ)
        res.append(duplicate)
    return res

def Gaussian2D((xs, ys), amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi):
    x, y = xs - x0, ys - y0
    xP = cosPhi * x - sinPhi * y
    yP = sinPhi * x + cosPhi * y
    res = amp * np.exp(-((xP) ** 2.0) / (2.0 * sigmaX ** 2.0) - \
                            ((yP) ** 2.0) / (2.0 * sigmaY ** 2.0)) + offset
    return res

def splevFast(x, tck):
    t, c, k = tck

    # if ext not in (0, 1, 2, 3):
    #    raise ValueError("ext = %s not in (0, 1, 2, 3) " % ext)

    # x = np.asarray(x)
    x = x.ravel()
    y, _ = interpolate._fitpack._spl_(x, 0, t, c, k, 0)

    # if ier == 10:
    #    raise ValueError("Invalid input data")
    # if ier == 1:
    #    raise ValueError("Found x value not in the domain")
    # if ier:
    #    raise TypeError("An error occurred")
    return y

def FitLeastsqScipy(errfunc, p0, **kwargs):
    method = kwargs.pop("method", "lm")
    if method != "lm":
        raise ValueError("Only method lm-supported")
    
    pfit, pcov, _, __, ___ = optimize.leastsq(errfunc, p0, full_output = 1, **kwargs)
    residual = errfunc(pfit)

    if (len(residual) > len(p0)) and pcov is not None:
        s_sq = (residual ** 2.0).sum() / (len(residual) - len(p0))
        pcov = pcov * s_sq
    else:
        pcov = np.inf

    error = [] 
    for i in range(len(pfit)):
        try:
            error.append(np.absolute(pcov[i][i]) ** 0.5)
        except:
            error.append(0.0)
    pfit_leastsq = pfit
    perr_leastsq = np.array(error) 
    return pfit_leastsq, perr_leastsq, residual
