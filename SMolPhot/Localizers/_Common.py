"""Implements general features required by localizers.

"""

import numpy as np
from SMolPhot.Components.BaseClasses import ParamsBaseClass

__all__ = []

                
class SimpleRegion():
    def __init__(self, offset, bbox):
        self.offset = offset
        self.bbox = bbox
                
#===============================================================================
# AddLocsWoDuplicates
#===============================================================================

def AddLocsWoDuplicates(locs, locsToAdd, minSameIterationDistXY, minDistXY, logger):
        # Modifies locs
        # Assumes, that locs contain no duplicates
        distFunc = lambda loc1, loc2: np.sqrt((loc1.x - loc2.x) ** 2.0 + (loc1.y - loc2.y) ** 2.0)
    
        duplicate = [False] * len(locsToAdd)        
        for i in range(len(locsToAdd) -1, -1, -1):
            # Check if locsToAdd[i] is duplicate
              
            # Check against old locations
            for j in range(len(locs)):
                if duplicate[i]:
                    # if already duplicate
                    break
                
                distXY = distFunc(locsToAdd[i], locs[j])
                if distXY < minDistXY:
                    logger.debug("Found duplicate location (cmp with old locs), removing.")
                    duplicate[i] = True
                    
            # Check against new locations
            for j in range(i):
                if duplicate[i]:
                    # if already duplicate
                    break
                
                distXY = distFunc(locsToAdd[i], locsToAdd[j])
                if distXY < minSameIterationDistXY:
                    logger.debug("Found duplicate location (cmp with new locs), removing.")
                    duplicate[i] = True
        
        locsToAddNew = [loc for i, loc in enumerate(locsToAdd) if not duplicate[i]]
        locs += locsToAddNew
        
        return locs, locsToAddNew
        

#===============================================================================
# CalcMinimumDistance
#===============================================================================

def CalcMinimumDistance(locs):
    
    for i in range(len(locs)):
        # find minimum distance
        minDistXY = np.inf
        for j in range(len(locs)):
            if i == j:
                continue
            distXY = np.sqrt((locs[i].x - locs[j].x) ** 2.0 + (locs[i].y - locs[j].y) ** 2.0)
            if distXY < minDistXY:
                minDistXY = distXY
        # save minimum distance
        minDistXY = min(10e-6, minDistXY) # TODO: 
        locs[i].minFitDistXY = minDistXY
    
    return locs

#===============================================================================
# LocalizerBase
#===============================================================================

class LocalizerBase(ParamsBaseClass):
    """Baseclass for all localizers.
    
    """
    
    def OptimizeForSeries(self, fseries):
        """Optimizes localizer parameters for given frame series.
        
        Args:
            fseries (FrameSeries): frame series
        
        """
        
        pass
    
    def FindMolecules(self, frame, psf, axialCalibrator):
        """Finds molecules in single frame.
        
        Args:
            frame (Frame): frame to analyze
            psf (PSF): PSF to use
        
        Returns:
            list of MoleculeLoc 
        
        """
        raise NotImplementedError()
    
    @property
    def name(self):
        raise NotImplementedError()

if __name__ == '__main__':
    pass