import logging
import numpy as np
from scipy import spatial
from _Common import PostprocessorBase
from SMolPhot.Components.BaseClasses import Param

__all__ = ["NearestMoleculePostprocessor"]

class NearestMoleculePostprocessor(PostprocessorBase):
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Postprocessor.NearestMoleculePostprocessor")
        paramsThis = [Param("_maxIterations", 10,
                        friendlyName = "Max iterations",
                        guiStep = 1,
                        guiLimits = (1, 1000)),
                      
                      Param("_nearestMolMaxDistXY", 60e-9,
                        friendlyName = "Nearest max dist XY",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 1e-9),
                      
                      Param("_minPercentOfMoleculesXY", 1.0,
                        friendlyName = "Min % of molecules (XY)",
                        guiStep = 0.05,
                        guiLimits = (0.0, 100)),
                      
                      Param("_nearestMolMaxDistXYZ", 60e-9,
                        friendlyName = "Nearest max dist XYZ",
                        guiSiPrefix = True,
                        guiSuffix = "m",
                        guiStep = 1e-9),
                      
                      Param("_minPercentOfMoleculesXYZ", 1.0,
                        friendlyName = "Min % of molecules (XYZ)",
                        guiStep = 0.05,
                        guiLimits = (0.0, 100.0)),
                      ]
        PostprocessorBase.__init__(self, paramsThis + params, **kwargs)
        
             
              
    def Apply(self, fseries, psf, groundTruth, locsToProcess):
        
        def CountNeighbours(coordinates, maxMolecules, maxDist):
            kdtree = spatial.cKDTree(coordinates)
            _, minDistIndices = kdtree.query(coordinates, k = maxMolecules, \
                                                    distance_upper_bound = maxDist)
            res = np.sum(minDistIndices < len(locsInput), axis = 1)
            return res
        
        #----------------------------------------------------------------------- 
        
        minNrOfMoleculesXY = max(1, int(round(0.01 * self._minPercentOfMoleculesXY * len(locsToProcess))))
        minNrOfMoleculesXYZ = max(1, int(round(0.01 * self._minPercentOfMoleculesXYZ * len(locsToProcess)))) 
        print minNrOfMoleculesXY, minNrOfMoleculesXYZ
        
        if not self._enabled or len(locsToProcess) <= 1:
            return locsToProcess
        self.logger.info("Apply to %d locs." % (len(locsToProcess)))
        
        locsInput = locsToProcess
        for iteration in range(self._maxIterations):
            self.logger.info("Iteration: %d" % (iteration))
            
            # Compose coordinates
            coordinatesXYZ = np.zeros((len(locsInput), 3), dtype = float)
            for i, loc in enumerate(locsInput):
                coordinatesXYZ[i, 0] = loc.x
                coordinatesXYZ[i, 1] = loc.y
                coordinatesXYZ[i, 2] = 0.0 if loc.z is None else loc.z
            coordinatesXY = coordinatesXYZ[:, [0, 1]]
            
            # Count neighbours
            countsXY = CountNeighbours(coordinatesXY, minNrOfMoleculesXY + 1, self._nearestMolMaxDistXY)
            countsXYZ = CountNeighbours(coordinatesXYZ, minNrOfMoleculesXYZ + 1, self._nearestMolMaxDistXYZ)
                
            # Discard
            locs = []
            for loc, countXY, countXYZ in zip(locsInput, countsXY, countsXYZ):
                if countXY < minNrOfMoleculesXY + 1:
                    self.logger.debug("Discarded, not enough XY points (%d)" % countXY)
                    continue
                
                if countXYZ < minNrOfMoleculesXYZ + 1:
                    self.logger.debug("Discarded, not enough XYZ points (%d)" % countXYZ)
                    continue
                locs.append(loc)
            
            self.logger.info("Discarded %d locs." % (len(locsInput) - len(locs)))    
            # Break if nothing happened
            if len(locs) == len(locsInput) or len(locs) <= 1:
                break
            
            locsInput = locs
        return locs

    @property
    def name(self):
        return "Nearest molecule postprocessor"
    
    