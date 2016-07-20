import math
import logging
import traceback
import numpy as np

from scipy import optimize
from collections import OrderedDict
from SMolPhot.Components.BaseClasses import Param
from _Common import PsfBase, MoleculeLoc, AxialCalibParam, FitLeastsqScipy, \
    Gaussian2D, Gaussian2DSym

# Try to import optional libraries
try:
    import GaussianFitter
except:
    print "Importing GaussianFitter failed"
    traceback.print_exc()

__all__ = ["GaussianPsf"]

#===============================================================================
# GaussianPsf
#===============================================================================

class GaussianPsf(PsfBase):

    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.PSFs.GaussianPsf")
        
        # Goodness functions
        self.goodnessFunctions = OrderedDict([
           ("amp * nrOfPixels / lsError", lambda loc: (abs(loc.amp) * float(loc.nrOfPixels)) / loc.lsError),
           ("('amp', 'fitStdOffset') (1.0, -2.5)", lambda loc: loc.amp ** 1.0 * loc.fitStdOffset ** -2.5),
           ("('amp', 'offset', 'fitStdOffset') (2.0, -0.5, -1.5)", lambda loc: loc.amp ** 2.0 * abs(loc.offset) ** -0.5 * loc.fitStdOffset ** -1.5),
           ("('photons', 'fitStdAmp', 'fitStdY0') (0.5, -1.5, -2.5)", lambda loc: loc.photons ** 0.5 * loc.fitStdAmp ** -1.5 * loc.fitStdY0 ** -2.5),
           ("distXY", lambda l: 1.0 / l.distXY),
           ("distXYZ", lambda l: 1.0 / np.sqrt(l.distXY ** 2.0 + l.distZ ** 2.0))
           ])
        
        # Error fonctions
        self.errorFunctions = OrderedDict([
           ("photons^0.5 * stdx^2 * minFitDistXY^-1.5", lambda loc: loc.photons ** 0.5 * loc.fitStdX0 ** 2.0 * loc.minFitDistXY ** -1.5),
           ("('photons', 'fitStdAmp', 'fitStdOffset') (1.5, -1.5, -0.5)", lambda loc: loc.photons ** -1.5 * loc.fitStdAmp ** 1.5 * loc.fitStdOffset ** 0.5)
           ])
        
        # Fitting libraries
        self._fittingLibrarys = OrderedDict([("SciPy leastsq", "scipyLeastsq"),
                                             ("SciPy curve_fit", "scipyCurveFit"),
                                             ("Fortran gaussian fitter", "GaussianFitter")])
        
        # Fitting methods
        self._fittingMethods = OrderedDict([("Levenberg-Marquardt", "lm"),
                                            (" Trust Region Reflective", "trf")])
        
        goodnessFuncsNames = self.goodnessFunctions.keys()
        errorFuncsNames = self.errorFunctions.keys()
        
        paramsThis = [
            # Symmetric
            
            Param("_symmetric", False,
                  friendlyName = "Symmetric (2D/3D)"),
            
            # Axial calibration
            Param("_sIntSigmaX", 0.0,
                  friendlyName = "Smooth sigmaX",
                  guiStep = 0.1,
                  guiLimits = (0.0, 1e12)),
            
            Param("_sIntSigmaY", 0.0,
                  friendlyName = "Smooth sigmaY",
                  guiStep = 0.1,
                  guiLimits = (0.0, 1e12)),
                      
            # Wobble
                       
            Param("_doWobbleCorrection", True,
                  friendlyName = "Wobble correction"),
            
            Param("_sIntDx", 0.0,
                  friendlyName = "Smooth wobble x",
                  guiStep = 0.1,
                  guiLimits = (0.0, 1e12)),
            
            Param("_sIntDy", 0.0,
                  friendlyName = "Smooth wobble y",
                  guiStep = 0.1,
                  guiLimits = (0.0, 1e12)),
                      
            # Fitting        
            
            Param("_fittingLibrary", "scipyLeastsq",
                  friendlyName = "Fitting library",
                  guiType = "list",
                  guiValues = self._fittingLibrarys),
            
            Param("_fittingMethod", "lm",
                  friendlyName = "Fitting method",
                  guiType = "list",
                  guiValues = self._fittingMethods),
                      
            Param("_fitXTol", 1.4901161193847656e-08,
                  friendlyName = "Fit x-tol"),

            Param("_fitFTol", 1.4901161193847656e-08,
                  friendlyName = "Fit f-tol"),

            Param("_fitGTol", 1.4901161193847656e-08,
                  friendlyName = "Fit g-tol"),
                      
            Param("_fitScale", False,
                  friendlyName = "Fit scale"),

            Param("_initialSigma", 300e-9,
                  friendlyName = "Initial sigma",
                  guiSiPrefix = True,
                  guiSuffix = "m",
                  guiStep = 10e-9),
            
            Param("_initialOffset", 0.0,
                  friendlyName = "Initial offset"),
            
            Param("_initialZ0", 0.0,
                  friendlyName = "Initial z",
                  guiSiPrefix = True,
                  guiSuffix = "m",
                  guiStep = 10e-9),
            
            Param("_phi", 0.0,
                  friendlyName = "Phi",
                  guiStep = 0.0001,
                  guiLimits = (-3.15, 3.15)),
            
            Param("_fitAreaShiftMaxIterations", 0,
                  friendlyName = "Max shift iterations",
                  guiStep = 1,
                  guiLimits = (0, 5)),
            
            # Goodness/Error functions
            
            Param("_goodnessFunc", goodnessFuncsNames[0],
                  friendlyName = "Goodness func",
                  guiType = "list",
                  guiValues = goodnessFuncsNames),
            
            Param("_errorFunc", errorFuncsNames[0],
                  friendlyName = "Error func",
                  guiType = "list",
                  guiValues = errorFuncsNames),
            ]
        
        PsfBase.__init__(self, params + paramsThis, **kwargs)
        self._Update()

    def SetParams(self, default = False, **kwargs):
        PsfBase.SetParams(self, default = default, **kwargs)
        if self._initDone:
            self._Update()

    def Calc(self, (xs, ys), amp, x0, y0, sigmaX, sigmaY, offset):
        return Gaussian2D((xs, ys), amp, x0, y0, sigmaX, sigmaY, offset, self.cosPhi, self.sinPhi)
        
    def CalcWithoutOffset(self, (xs, ys), amp, x0, y0, sigmaX, sigmaY, _):
        return Gaussian2D((xs, ys), amp, x0, y0, sigmaX, sigmaY, 0.0, self.cosPhi, self.sinPhi)

    def Fit(self, frame, initialCoord, fitPixels, fitMode = "z", initialMaxDist = np.inf, \
            initialGuess = None, fitShiftPixels = (0, 0), fitShiftIteration = 0):
                
        # Get data around initialCoord
        initialPixels = frame.GetCoordsPixels(initialCoord)
        fitCenterPixels = (initialPixels[0] + fitShiftPixels[0], initialPixels[1] + fitShiftPixels[1])
        if fitShiftIteration == 0:
            fitCenterPixels = initialPixels
        else:
            fitCenterPixels = frame.CropCoord((initialPixels[0] + fitShiftPixels[0], \
                                               initialPixels[1] + fitShiftPixels[1]))
        _, coordsFlatten, dataFlatten = frame.GetAround(fitCenterPixels, fitPixels, flatten = True)
        iAmp = frame._data[fitCenterPixels]

        # Do fitting
        try:
            initial = initialGuess
            if fitMode == "sigma":
                fitParamValues, fitStds, residual = self._DoFittingSigma(iAmp, initialCoord, \
                    coordsFlatten, dataFlatten, initial, **self.fitExtraParams)                
            elif fitMode == "z":
                fitParamValues, fitStds, residual = self._DoFittingZ(iAmp, initialCoord, \
                    coordsFlatten, dataFlatten, initial, fitPixels, **self.fitExtraParams)     
            else:
                raise NotImplementedError()
            
        except RuntimeError:
            traceback.print_exc()
            return None
            
        # Extract params
        amp, x0, y0, z0, sigmaX, sigmaY, offset = fitParamValues
        stdAmp, stdX0, stdY0, stdZ0, stdSigmaX, stdSigmaY, stdOffset = fitStds
            
        # Check if is close to the initial coord
        if np.sqrt((x0 - initialCoord[0]) ** 2.0 + (y0 - initialCoord[1]) ** 2.0) > initialMaxDist:
            return None
        
        # Check if gaussian maximum is in the same pixel as expected
        fittedMaximumPixels = frame.GetCoordsPixels((x0, y0))
        if fitShiftIteration < self._fitAreaShiftMaxIterations and \
            fittedMaximumPixels != fitCenterPixels:
            
            self.logger.debug("Fit maximum in wrong pixel: %s %s" % \
                              (fittedMaximumPixels, fitCenterPixels))
            shiftPixels = (fittedMaximumPixels[0] - fitCenterPixels[0],
                           fittedMaximumPixels[1] - fitCenterPixels[1])
            # Initial
            if fitMode == "sigma":
                initial = (amp, x0, y0, sigmaX, sigmaY, offset)
            elif fitMode == "z":
                initial = (amp, x0, y0, z0, offset)
            else:
                raise NotImplementedError()
            
            # Do new fitting from shifted position
            locFitNew = self.Fit(frame, initialCoord, fitPixels, \
                            fitMode = fitMode, \
                            initialMaxDist = initialMaxDist, \
                            initialGuess = initial, \
                            fitShiftPixels = shiftPixels, \
                            fitShiftIteration = fitShiftIteration + 1)
            
            # If found better
            if locFitNew is not None:
                return locFitNew
            
        # Wobble correction
        xCor, yCor = x0, y0
        if self._doWobbleCorrection and fitMode == "z" and z0 is not None:
            dx = self._calibDx.GetInterpolatedValue(z0, fitPixels)
            dy = self._calibDy.GetInterpolatedValue(z0, fitPixels)
            self.logger.debug("Wobble correction: dx %.3f nm, dy %.3f nm at z %.3f nm" % \
                              (1e9 * dx, 1e9 * dy, 1e9 * z0))
            xCor, yCor = x0 + dx[0], y0 + dy[0]

        # Result
        photons = 2.0 * np.pi * amp * sigmaX * sigmaY / (frame.pixelSize ** 2.0)
        res = MoleculeLoc(coord = (xCor, yCor, z0),
                          coordPixels = frame.GetCoordsPixels((xCor, yCor)),
                          frameNr = frame.nr,
                          initial = initialCoord,
                          bestFitParams = (amp, x0, y0, sigmaX, sigmaY, offset),
                          residual = residual,
                          fitStds = fitStds,
                          photons = photons,
                          amp = amp,
                          sigmaX = sigmaX,
                          sigmaY = sigmaY,
                          offset = offset,
                          fitStdAmp = stdAmp,
                          fitStdX0 = stdX0,
                          fitStdY0 = stdY0,
                          fitStdZ0 = stdZ0,
                          fitStdSigmaX = stdSigmaX,
                          fitStdSigmaY = stdSigmaY,
                          fitStdOffset = stdOffset,
                          fitPixels = fitPixels)
        
        # Calculate goondess end error values
        res.goodnessValue = self.goodnessFunc(res)
        res.errorValue = self.errorFunc(res)
        return res
    
    def SetAxialCalibPoints(self, points, useStds = False):
        self.logger.info("SetAxialCalibPoints %d" % (len(points)))
        
        if len(points) < 1:
            self._calibSigmaX, self._calibSigmaY = None, None
            self._calibDx, self._calibDy = None, None
            self._hasCalibrationData = False
            return
            
        self._calibSigmaX = AxialCalibParam(points, lambda loc: loc.sigmaX, \
                                            self._sIntSigmaX, useStds = useStds)
        
        self._calibSigmaY = AxialCalibParam(points, lambda loc: loc.sigmaY, \
                                            self._sIntSigmaY, useStds = useStds)
        
        self._calibDx = AxialCalibParam(points, lambda loc: loc.startPositionOffset[0],
                                        self._sIntDx, useStds = useStds)
        
        self._calibDy = AxialCalibParam(points, lambda loc: loc.startPositionOffset[1],
                                        self._sIntDy, useStds = useStds)
        
        self._hasCalibrationData = True
        
    def UpdateCalibInterpolation(self):
        if not self.hasCalibrationData:
            return
        
        self._calibSigmaX.SetSmoothing(self._sIntSigmaX)
        self._calibSigmaY.SetSmoothing(self._sIntSigmaY)
        self._calibDx.SetSmoothing(self._sIntDx)
        self._calibDy.SetSmoothing(self._sIntDy)

    @property
    def wobbleCalibParamValues(self):
        return [("dx", self._calibDx), ("dy", self._calibDy)]
    
    @property
    def zCalibParamValues(self):
        return [("sigmaX", self._calibSigmaX), ("sigmaY", self._calibSigmaY)]
    
    @property
    def goodnessFunc(self):
        return self.goodnessFunctions[self._goodnessFunc]
    
    @property
    def errorFunc(self):
        return self.errorFunctions[self._errorFunc]

    @property
    def name(self):
        return "Gaussian"
    
    def _Update(self):
        self.sinPhi, self.cosPhi = math.sin(self._phi), math.cos(self._phi)
        self.fitExtraParams = {"method": self._fittingMethod, "xtol": self._fitXTol, \
                               "ftol": self._fitFTol, "gtol": self._fitGTol}
        self.UpdateCalibInterpolation()
    
    def _DoFittingSigma(self, iAmp, initialCoord, coords, data, initial, **kwargs):
        # Initial
        initial = [iAmp, 
                   initialCoord[0], 
                   initialCoord[1], 
                   self._initialSigma, 
                   self._initialSigma, 
                   self._initialOffset] if initial is None else initial
        
        # Scale
        if self._fitScale:
            self.fitExtraParams["x_scale"] = [initial[0], 3e-6, 3e-6, 300e-9, 300e-9, 10.0]
            
        # Differnet libraries
        if self._fittingLibrary == "scipyLeastsq":
            if self._symmetric:
                fitRes = Fit2DGaussianSigmaLeastSq(coords, data, initial, **kwargs)
            else:
                fitRes = Fit3DGaussianSigmaLeastSq(coords, data, initial, self.cosPhi, self.sinPhi, **kwargs)
            fitParamValues, fitStds, residual = fitRes
                
        elif self._fittingLibrary == "scipyCurveFit":
            if self._symmetric:
                raise NotImplementedError()
            else:
                fitRes = Fit3DGaussianSigmaCurveFit(coords, data, initial, self.cosPhi, self.sinPhi, **kwargs)
            fitParamValues, fitStds, residual = fitRes
        elif self._fittingLibrary == "GaussianFitter":
            method = kwargs.pop("method", "lm")
            if method != "lm":
                raise ValueError("Only Levenberg-Marquardt is supporeted by this library")
            if self._symmetric:
                initialNew = RemoveSigma(initial)
                fitRes = GaussianFitter.Fit2dGaussSigmaSym(coords[0], coords[1], data, initialNew, **kwargs)
                fitParamValues, pcov, residual = fitRes
                fitParamValues = DuplicateSigma(fitParamValues)
                fitStds = DuplicateSigma(np.sqrt(np.diag(pcov)))
            else:
                fitRes = GaussianFitter.Fit2dGaussSigma(coords[0], coords[1], data, initial, phi = self._phi, **kwargs)
                fitParamValues, pcov, residual = fitRes
                fitStds = np.sqrt(np.diag(pcov))
        else:
            raise NotImplementedError()
        
        # Extract params
        amp, x0, y0, sigmaX, sigmaY, offset = fitParamValues
        stdAmp, stdX0, stdY0, stdSigmaX, stdSigmaY, stdOffset = fitStds
            
        # Return params
        return (amp, x0, y0, None, sigmaX, sigmaY, offset), \
            (stdAmp, stdX0, stdY0, None, stdSigmaX, stdSigmaY, stdOffset), \
            residual

    def _DoFittingZ(self, iAmp, initialCoord, coords, data, initial, fitPixels, **kwargs):
        # Initial
        initial = [iAmp, 
                   initialCoord[0], 
                   initialCoord[1], 
                   self._initialZ0, 
                   self._initialOffset] if initial is None else initial

        # Scale  
        if self._fitScale:
            self.fitExtraParams["x_scale"] = [initial[0], 3e-6, 3e-6, 1e-6, 10.0]
                
        sigmaXFunc = lambda z: self._calibSigmaX.GetInterpolatedValue(z, fitPixels)
        sigmaYFunc = lambda z: self._calibSigmaY.GetInterpolatedValue(z, fitPixels)
        
        if self._fittingLibrary == "scipyLeastsq":
            if self._symmetric:
                fitRes = Fit2DGaussianZLeastSq(coords, data, initial, sigmaXFunc, **kwargs)
            else:
                fitRes = Fit3DGaussianZLeastSq(coords, data, initial, self.cosPhi, \
                    self.sinPhi, sigmaXFunc, sigmaYFunc, **kwargs)
            fitParamValues, fitStds, residual = fitRes
        elif self._fittingLibrary == "scipyCurveFit":
            if self._symmetric:
                raise NotImplementedError()
            else:
                fitRes = Fit3DGaussianZCurveFit(coords, data, initial, self.cosPhi, \
                    self.sinPhi, sigmaXFunc, sigmaYFunc, **kwargs)
            fitParamValues, fitStds, residual = fitRes
        elif self._fittingLibrary == "GaussianFitter":
            method = kwargs.pop("method", "lm")
            if method != "lm":
                raise ValueError("Only Levenberg-Marquardt is supporeted by this library")
            tckSigmaX = self._calibSigmaX.GetData(fitPixels, error = False)["interpCoefs"]
            tckSigmaY = self._calibSigmaY.GetData(fitPixels, error = False)["interpCoefs"]
            if self._symmetric:
                fitRes = GaussianFitter.Fit2dGaussZSym(coords[0], coords[1], data, \
                    initial, tckSigmaX, **kwargs)
            else:
                fitRes = GaussianFitter.Fit2dGaussZ(coords[0], coords[1], data, \
                    initial, tckSigmaX, tckSigmaY, phi = self._phi, **kwargs)
            fitParamValues, pcov, residual = fitRes
            fitStds = np.sqrt(np.diag(pcov))
        else:
            raise NotImplementedError()
        
        # Extract params
        amp, x0, y0, z0, offset = fitParamValues
        stdAmp, stdX0, stdY0, stdZ0, stdOffset = fitStds
            
        # Return params
        return (amp, x0, y0, z0, sigmaXFunc(z0), sigmaYFunc(z0), offset), \
            (stdAmp, stdX0, stdY0, stdZ0, None, None, stdOffset), \
            residual
            
    
#===============================================================================
# Methods
#===============================================================================
# Scipy.leastsq 2D

def RemoveSigma(array):
    return (array[0], array[1], array[2], array[3], array[5])

def DuplicateSigma(array):
    return (array[0], array[1], array[2], array[3], array[3], array[4])

def Fit2DGaussianSigmaLeastSq(coords, data, initial, **kwargs):        
    def FitSigmaFunc(params):
        amp, x0, y0, sigma, offset = params
        calculated = Gaussian2DSym(coords, amp, x0, y0, sigma, offset)
        r = data - calculated
        return r
    
    initialNew = RemoveSigma(initial)
    fitParamValues, fitStds, residual = FitLeastsqScipy(FitSigmaFunc, initialNew, **kwargs)
    return DuplicateSigma(fitParamValues), DuplicateSigma(fitStds), residual

def Fit2DGaussianZLeastSq(coords, data, initial, sigmaFunc, **kwargs):   
    def FitZFunc(params):
        amp, x0, y0, z0, offset = params
        sigma = sigmaFunc(z0)
        calculated = Gaussian2DSym(coords, amp, x0, y0, sigma, offset)
        r = data - calculated
        return r
    
    fitParamValues, fitStds, residual = FitLeastsqScipy(FitZFunc, initial, **kwargs)
    return fitParamValues, fitStds, residual

def Fit3DGaussianSigmaCurveFit(coords, data, initial, cosPhi, sinPhi, **kwargs):        
    def Func(coords, amp, x0, y0, sigmaX, sigmaY, offset):
        calculated = Gaussian2D(coords, amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi)
        return calculated
    
    if kwargs["method"] == "lm":
        kwargs["diag"] = kwargs.pop("x_scale", None)
    
    fitParamValues, pcov = optimize.curve_fit(Func, coords, data, initial, **kwargs)
    fitStds = np.sqrt(np.diag(pcov))
    residual = data - Func(coords, *fitParamValues)
    return fitParamValues, fitStds, residual

def Fit3DGaussianZCurveFit(coords, data, initial, cosPhi, sinPhi, sigmaXFunc, sigmaYFunc, **kwargs):        
    def Func(coords, amp, x0, y0, z0, offset):
        sigmaX = sigmaXFunc(z0)
        sigmaY = sigmaYFunc(z0)
        calculated = Gaussian2D(coords, amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi)
        return calculated
    
    if kwargs["method"] == "lm":
        kwargs["diag"] = kwargs.pop("x_scale", None)
        
    fitParamValues, pcov = optimize.curve_fit(Func, coords, data, initial, **kwargs)
    fitStds = np.sqrt(np.diag(pcov))
    residual = data - Func(coords, *fitParamValues)
    return fitParamValues, fitStds, residual


# Scipy.leastsq 3D

def Fit3DGaussianSigmaLeastSq(coords, data, initial, cosPhi, sinPhi, **kwargs):        
    def FitSigmaFunc(params):
        amp, x0, y0, sigmaX, sigmaY, offset = params
        calculated = Gaussian2D(coords, amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi)
        r = data - calculated
        return r
    
    fitParamValues, fitStds, residual = FitLeastsqScipy(FitSigmaFunc, initial, **kwargs)
    return fitParamValues, fitStds, residual

def Fit3DGaussianZLeastSq(coords, data, initial, cosPhi, sinPhi, sigmaXFunc, sigmaYFunc, **kwargs):   
    def FitZFunc(params):
        amp, x0, y0, z0, offset = params
        sigmaX = sigmaXFunc(z0)
        sigmaY = sigmaYFunc(z0)
        calculated = Gaussian2D(coords, amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi)
        r = data - calculated
        return r
    
    fitParamValues, fitStds, residual = FitLeastsqScipy(FitZFunc, initial, **kwargs)
    return fitParamValues, fitStds, residual
      
if __name__ == "__main__":
    pass
