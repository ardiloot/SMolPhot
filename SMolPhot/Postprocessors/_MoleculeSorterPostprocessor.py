import logging
import numpy as np

from scipy import spatial
from _Common import PostprocessorBase
from SMolPhot.Components.BaseClasses import Param

__all__ = ["MoleculeSorterPostprocessor"]

class MoleculeSorterPostprocessor(PostprocessorBase):
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Postprocessor.MoleculeSorterPostprocessor")
        paramsThis = [Param("_percentageToDiscard", 5.0, \
                        friendlyName = "To discard", \
                        guiSuffix = " %", \
                        guiStep = 1.0, \
                        guiLimits = (0, 100)),
                      
                      Param("_nearestMolMaxDistXY",
                        {"checked": False, "value": 60e-9},
                        friendlyName = "Nearest max dist XY",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 1e-9),
                      
                      Param("_nearestMolMaxDistXYZ",
                        {"checked": False, "value": 60e-9},
                        friendlyName = "Nearest max dist XYZ",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 1e-9),
                      
                      Param("_errorFuncPercentage",
                        {"checked": False, "value": 5.0},
                        friendlyName = "Error func threshold (%)",
                        guiLimits = (0.0, 100.0),
                        guiStep = 0.4),
                      
                      Param("_maxFitStdX0",
                        {"checked": False, "value": 200e-9},
                        friendlyName = "Max std x",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 1e-9),
                      
                      Param("_maxFitStdY0",
                        {"checked": False, "value": 200e-9},
                        friendlyName = "Max std y",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 1e-9),
                      
                      Param("_maxFitStdZ0",
                        {"checked": False, "value": 200e-9},
                        friendlyName = "Max std z",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 1e-9),
                      
                      Param("_maxFitStdAmp",
                        {"checked": False, "value": 10},
                        friendlyName = "Max std amp",
                        guiStep = 0.1),
                      
                      Param("_maxFitStdOffset",
                        {"checked": False, "value": 10},
                        friendlyName = "Max std offset",
                        guiStep = 0.1),
                      
                      Param("_minAmp",
                        {"checked": False, "value": 20},
                        friendlyName = "Min amp",
                        guiStep = 0.1),
                      
                      Param("_minPhotons",
                        {"checked": False, "value": 20},
                        friendlyName = "Min photons",
                        guiStep = 0.1),
                      
                      Param("_maxZ",
                        {"checked": False, "value": 1e-6},
                        friendlyName = "Max z",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 10e-9),
                      
                      Param("_minZ",
                        {"checked": False, "value": 1e-6},
                        friendlyName = "Min z",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 10e-9),
                      
                      Param("_maxOffset",
                        {"checked": False, "value": 100},
                        friendlyName = "Max offset",
                        guiStep = 0.1),
                      
                      Param("_minOffset",
                        {"checked": False, "value": -100},
                        friendlyName = "Min offset",
                        guiStep = 0.1),
                                            
                      Param("_minFitDistXY",
                        {"checked": False, "value": 100e-9},
                        friendlyName = "Min fit dist",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 10e-9)
                      ]
        PostprocessorBase.__init__(self, paramsThis + params, **kwargs)
        
    def _UpdateNeighboursDists(self, locs):
        
        def GetNeighboursDistances(coordinates):
            kdtree = spatial.cKDTree(coordinates)
            minDists, minDistIndices = kdtree.query(coordinates, k = 2)
            if np.any(abs(minDists[:, 0]) > 1e-16):
                raise RuntimeError("Closest-neighbour should always be at zero distance (itself).") 
            
            # Discard closest (itself)
            minDists = minDists[:, 1]
            minDistIndices = minDistIndices[:, 1]
            return minDists, minDistIndices
        
        if len(locs) < 2:
            for loc in locs:
                loc.finalNearestMolDistXY = np.inf
                loc.finalNearestMolDistXYZ = np.inf
            return
        
        coordinatesXYZ = np.zeros((len(locs), 3), dtype = float)
        for i, loc in enumerate(locs):
            coordinatesXYZ[i, 0] = loc.x
            coordinatesXYZ[i, 1] = loc.y
            coordinatesXYZ[i, 2] = 0.0 if loc.z is None else loc.z
        coordinatesXY = coordinatesXYZ[:, [0, 1]]
      
        minDistXY, _ = GetNeighboursDistances(coordinatesXY)
        minDistXYZ, _ = GetNeighboursDistances(coordinatesXYZ)
        
        for i, loc in enumerate(locs):
            loc.finalNearestMolDistXY = minDistXY[i]
            loc.finalNearestMolDistXYZ = minDistXYZ[i]
             
    def _UpdatedGoodnessAndErrorValues(self, psf, locs):
        for loc in locs:
            loc.goodnessValue = psf.goodnessFunc(loc)
            loc.errorValue = psf.errorFunc(loc)
              
    def Apply(self, fseries, psf, groundTruth, locsToProcess):
        if not self._enabled or len(locsToProcess) <= 0:
            return locsToProcess
        # Update locs precalculated values
        self._UpdatedGoodnessAndErrorValues(psf, locsToProcess)
        self._UpdateNeighboursDists(locsToProcess)
        
        # Build limitations
        limitFuncs = []
        
        if self._nearestMolMaxDistXY["checked"]:
            limitFuncs += [lambda loc: loc.finalNearestMolDistXY < self._nearestMolMaxDistXY["value"]]
        
        if self._nearestMolMaxDistXYZ["checked"]:
            limitFuncs += [lambda loc: loc.finalNearestMolDistXYZ < self._nearestMolMaxDistXYZ["value"]]
        
        if self._maxFitStdX0["checked"]:
            limitFuncs += [lambda loc: loc.fitStdX0 < self._maxFitStdX0["value"]]

        if self._maxFitStdY0["checked"]:
            limitFuncs += [lambda loc: loc.fitStdY0 < self._maxFitStdY0["value"]]
            
        if self._maxFitStdZ0["checked"]:
            limitFuncs += [lambda loc: loc.fitStdZ0 < self._maxFitStdZ0["value"]]
        
        if self._maxFitStdAmp["checked"]:
            limitFuncs += [lambda loc: loc.fitStdAmp < self._maxFitStdAmp["value"]]
        
        if self._maxFitStdOffset["checked"]:
            limitFuncs += [lambda loc: loc.fitStdOffset < self._maxFitStdOffset["value"]]
       
        if self._minAmp["checked"]:
            limitFuncs += [lambda loc: loc.amp > self._minAmp["value"]]
       
        if self._minPhotons["checked"]:
            limitFuncs += [lambda loc: loc.photons > self._minPhotons["value"]]
            
        if self._maxZ["checked"]:
            limitFuncs += [lambda loc: loc.z < self._maxZ["value"]]
       
        if self._minZ["checked"]:
            limitFuncs += [lambda loc: loc.z > self._minZ["value"]]
            
        if self._maxOffset["checked"]:
            limitFuncs += [lambda loc: loc.offset < self._maxOffset["value"]]
       
        if self._minOffset["checked"]:
            limitFuncs += [lambda loc: loc.offset > self._minOffset["value"]]
            
        if self._minFitDistXY["checked"]:
            limitFuncs += [lambda loc: loc.minFitDistXY > self._minFitDistXY["value"]]
            
        # Selects bad points based on psf error function
        if self._errorFuncPercentage["checked"] and len(locsToProcess) > 0:
            errorValues = sorted([loc.errorValue for loc in locsToProcess], key = lambda k: -k)
            errorThreshold = errorValues[int(np.around(0.01 * self._errorFuncPercentage["value"] * (len(errorValues) - 1)))] 
            limitFuncs += [lambda loc: psf.errorFunc(loc) < errorThreshold]
            
        # Build total-goodness function
        goodnessFunc = lambda loc: [f(loc) for f in limitFuncs] + [loc.goodnessValue]
        locs = sorted(locsToProcess, key = goodnessFunc, reverse = True)

        # Discard less good points
        lastMolNr = max(1, np.around(0.01 * (100.0 - self._percentageToDiscard) * len(locs)).astype(int))
        resLocs = locs[:lastMolNr]        
        return resLocs

    @property
    def name(self):
        return "Molecule sorter postprocessor"
    
    