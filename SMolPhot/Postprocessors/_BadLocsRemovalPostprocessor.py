import logging
import numpy as np

from collections import OrderedDict
from _Common import PostprocessorBase
from SMolPhot.Components.BaseClasses import Param


__all__ = ["BadLocsRemovalPostprocessor"]

class BadLocsRemovalPostprocessor(PostprocessorBase):
    
    @staticmethod
    def GetBadnessFunctions():
        # Temporary, separate function good for debugging
        errorFunctions = OrderedDict()
        errorFunctions["MT0.N1.HD-3D-AS_15"] = [
            (0.003333, lambda loc: loc.photons ** -3.0 * loc.absoffset ** 1.0 * loc.fitStdY0 ** 1.0),
            (0.006667, lambda loc: loc.sigmaX ** 2.0 * loc.minFitDistXY ** -3.0 * loc.nrOfPixels ** -3.0),
            (0.010000, lambda loc: loc.lsError ** 1.0 * loc.photons ** 2.0 * loc.fitStdX0 ** 3.0),
            (0.013333, lambda loc: loc.photons ** -1.0 * loc.minFitDistXY ** -2.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.016667, lambda loc: loc.sigmaY ** 3.0 * loc.minFitDistXY ** -3.0 * loc.nrOfPixels ** -2.0),
            (0.020000, lambda loc: loc.lsError ** -1.0 * loc.sigmaX ** -3.0 * loc.fitStdX0 ** 1.0),
            (0.023333, lambda loc: loc.lsError ** 1.0 * loc.absz ** 3.0 * loc.sigmaY ** -3.0),
            (0.026667, lambda loc: loc.fitStdZ0 ** 1.0 * loc.minFitDistXY ** -3.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.030000, lambda loc: loc.absoffset ** 3.0 * loc.fitStdAmp ** -3.0 * loc.fitStdZ0 ** 2.0),
            (0.033333, lambda loc: loc.sigmaY ** 1.0 * loc.minFitDistXY ** -3.0 * loc.distFitWeightedCentroid ** 1.0),
            (0.036667, lambda loc: loc.sigmaX ** -1.0 * loc.fitStdAmp ** -2.0 * loc.fitStdOffset ** 2.0),
            (0.040000, lambda loc: loc.lsError ** -2.0 * loc.fitStdY0 ** 1.0 * loc.minFitDistXY ** -2.0),
            (0.043333, lambda loc: loc.sigmaY ** -2.0 * loc.fitStdY0 ** 1.0 * loc.minFitDistXY ** -1.0),
            (0.046667, lambda loc: loc.sigmaX ** -3.0 * loc.fitStdX0 ** 1.0 * loc.distFitWeightedCentroid ** 2.0),
            (0.050000, lambda loc: loc.sigmaX ** -3.0 * loc.fitStdX0 ** 3.0 * loc.fitStdZ0 ** -2.0),
            ]
        
        errorFunctions["MT0.N2.HD-3D-AS_15"] = [
            (0.003333, lambda loc: loc.absz ** 3.0 * loc.fitStdAmp ** 3.0 * loc.fitStdZ0 ** -2.0),
            (0.006667, lambda loc: loc.lsError ** 1.0 * loc.absz ** 3.0 * loc.sigmaY ** -3.0),
            (0.010000, lambda loc: loc.lsError ** -3.0 * loc.fitStdZ0 ** 2.0 * loc.minFitDistXY ** -2.0),
            (0.013333, lambda loc: loc.amp ** -3.0 * loc.sigmaY ** 1.0 * loc.fitStdZ0 ** 3.0),
            (0.016667, lambda loc: loc.lsError ** 1.0 * loc.absz ** 3.0 * loc.nrOfPixels ** 1.0),
            (0.020000, lambda loc: loc.fitStdX0 ** 2.0 * loc.minFitDistXY ** -3.0 * loc.nrOfPixels ** 1.0),
            (0.023333, lambda loc: loc.photons ** -2.0 * loc.minFitDistXY ** -3.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.026667, lambda loc: loc.absz ** 1.0 * loc.minFitDistXY ** -3.0 * loc.area ** -1.0),
            (0.030000, lambda loc: loc.lsError ** 1.0 * loc.absz ** 3.0 * loc.fitStdOffset ** 1.0),
            (0.033333, lambda loc: loc.absz ** -1.0 * loc.absoffset ** 2.0 * loc.fitStdOffset ** 3.0),
            (0.036667, lambda loc: loc.amp ** -2.0 * loc.minFitDistXY ** -3.0 * loc.area ** 1.0),
            (0.040000, lambda loc: loc.photons ** -3.0 * loc.fitStdX0 ** -1.0 * loc.fitStdY0 ** 3.0),
            (0.043333, lambda loc: loc.absz ** 2.0 * loc.fitStdAmp ** 1.0 * loc.distFitWeightedCentroid ** 1.0),
            (0.046667, lambda loc: loc.absoffset ** 3.0 * loc.fitStdAmp ** -2.0 * loc.fitStdZ0 ** 3.0),
            (0.050000, lambda loc: loc.sigmaY ** 2.0 * loc.nrOfPixels ** -3.0 * loc.area ** -1.0),
            ]
        
        errorFunctions["MT0.N1.LD-3D-AS_15"] = [
            (0.003333, lambda loc: loc.amp ** -2.0 * loc.absoffset ** 3.0 * loc.fitStdOffset ** -1.0),
            (0.006667, lambda loc: loc.absoffset ** 2.0 * loc.fitStdY0 ** 3.0 * loc.minFitDistXY ** -3.0),
            (0.010000, lambda loc: loc.photons ** -3.0 * loc.fitStdY0 ** 3.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.013333, lambda loc: loc.photons ** 1.0 * loc.amp ** -1.0 * loc.nrOfPixels ** -1.0),
            (0.016667, lambda loc: loc.photons ** -2.0 * loc.absoffset ** 3.0 * loc.fitStdY0 ** 3.0),
            (0.020000, lambda loc: loc.amp ** -1.0 * loc.minFitDistXY ** -1.0 * loc.distFitWeightedCentroid ** 2.0),
            (0.023333, lambda loc: loc.photons ** -3.0 * loc.absoffset ** 2.0 * loc.fitStdOffset ** 1.0),
            (0.026667, lambda loc: loc.photons ** 3.0 * loc.fitStdZ0 ** 2.0 * loc.nrOfPixels ** 1.0),
            (0.030000, lambda loc: loc.lsError ** -2.0 * loc.absz ** 3.0 * loc.sigmaY ** 2.0),
            (0.033333, lambda loc: loc.absoffset ** 3.0 * loc.fitStdX0 ** 2.0 * loc.fitStdY0 ** 2.0),
            (0.036667, lambda loc: loc.photons ** -2.0 * loc.absz ** 2.0 * loc.minFitDistXY ** -2.0),
            (0.040000, lambda loc: loc.amp ** -3.0 * loc.fitStdAmp ** 1.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.043333, lambda loc: loc.lsError ** 2.0 * loc.fitStdX0 ** 1.0 * loc.minFitDistXY ** 1.0),
            (0.046667, lambda loc: loc.amp ** -2.0 * loc.sigmaY ** 1.0 * loc.area ** 1.0),
            (0.050000, lambda loc: loc.sigmaX ** -2.0 * loc.fitStdY0 ** 2.0 * loc.minFitDistXY ** -1.0),
            ]
        
        errorFunctions["MT0.N2.LD-3D-AS_15"] = [
            (0.003333, lambda loc: loc.photons ** -2.0 * loc.absoffset ** 3.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.006667, lambda loc: loc.amp ** -3.0 * loc.minFitDistXY ** -1.0 * loc.area ** 3.0),
            (0.010000, lambda loc: loc.lsError ** 3.0 * loc.absz ** 3.0 * loc.fitStdAmp ** 1.0),
            (0.013333, lambda loc: loc.photons ** -3.0 * loc.absoffset ** 3.0 * loc.minFitDistXY ** -3.0),
            (0.016667, lambda loc: loc.sigmaY ** 2.0 * loc.distFitWeightedCentroid ** 1.0 * loc.area ** -2.0),
            (0.020000, lambda loc: loc.absoffset ** 3.0 * loc.fitStdY0 ** 2.0 * loc.fitStdOffset ** -2.0),
            (0.023333, lambda loc: loc.photons ** -3.0 * loc.absz ** 3.0 * loc.nrOfPixels ** -2.0),
            (0.026667, lambda loc: loc.amp ** -3.0 * loc.absoffset ** 1.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.030000, lambda loc: loc.absoffset ** 1.0 * loc.fitStdY0 ** 3.0 * loc.minFitDistXY ** -3.0),
            (0.033333, lambda loc: loc.photons ** -1.0 * loc.absoffset ** 3.0 * loc.fitStdX0 ** 3.0),
            (0.036667, lambda loc: loc.absoffset ** 1.0 * loc.fitStdZ0 ** 1.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.040000, lambda loc: loc.absz ** 1.0 * loc.fitStdX0 ** 1.0 * loc.fitStdOffset ** -2.0),
            (0.043333, lambda loc: loc.lsError ** 2.0 * loc.photons ** -2.0 * loc.distFitWeightedCentroid ** 3.0),
            (0.046667, lambda loc: loc.photons ** -3.0 * loc.minFitDistXY ** 1.0 * loc.distFitWeightedCentroid ** 1.0),
            (0.050000, lambda loc: loc.lsError ** 1.0 * loc.photons ** -2.0 * loc.absoffset ** 2.0),
            ]
        return errorFunctions
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Postprocessor.BadLocsRemovalPostprocessor")

        paramsThis = [Param("_badnessFunc", "MT0_N1_HD_AS_10",
                        friendlyName = "Error func",
                        guiType = "list",
                        guiValues = self.GetBadnessFunctions().keys()),
                        
                      Param("_removePercentage", 1.0,
                        friendlyName = "To remove (%)",
                        guiLimits = (0.0, 100.0),
                        guiStep = 0.1),
                      ]
        PostprocessorBase.__init__(self, paramsThis + params, **kwargs)
        
    def Apply(self, fseries, psf, groundTruth, locsToProcess):
        # Save input for later use
        self._locsToProcess = locsToProcess[:]
        
        if not self._enabled or len(locsToProcess) <= 1:
            return locsToProcess
        self.logger.info("Apply to %d locs." % (len(locsToProcess)))

        badnessFuncData = self.GetBadnessFunctions()[self._badnessFunc]
        totalPartToRemove = 0.01 * self._removePercentage
        
        locs = locsToProcess[:]
        lowerLimit = 0.0
        nrOfAlreadyRemoved = 0
        for i, (upperLimit, errorFunc) in enumerate(badnessFuncData):
            if totalPartToRemove <= lowerLimit:
                # Already removed enough
                break
            
            # Nr of points to remove
            partToRemove = upperLimit if totalPartToRemove >= upperLimit else totalPartToRemove
            if i + 1 == len(badnessFuncData):
                # Last error func, remove as many as neccesary
                partToRemove = totalPartToRemove
            nrOfPointsToRemove = int(round(partToRemove * float(len(locsToProcess)))) - nrOfAlreadyRemoved
            
            # Sort locs by errorFunc
            locs.sort(key = errorFunc, reverse = True)
            
            # Remove points
            pointsToRemoveSet = set(locs[:nrOfPointsToRemove])
            newLocs = [loc for loc in locs if loc not in pointsToRemoveSet] 
            
            # Update variables for next            
            lowerLimit = upperLimit
            nrOfAlreadyRemoved += nrOfPointsToRemove
            locs = newLocs

        
        self.logger.info("Apply done.")
        return locs
        
    @property
    def name(self):
        return "Bad molecule removal postprocessor"

if __name__ == "__main__":
    pass