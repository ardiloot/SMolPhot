"""This module implements classes for ground-truth information and for
verification of the results. 

"""

import logging
import numpy as np

from scipy import optimize, spatial
from collections import defaultdict
from SMolPhot.PSFs import GetLocsCoordsArray
from BaseClasses import ParamsBaseClass, Param

__all__ = ["GroundTruthSettings", "GroundTruthStats"]


def GetStatsStr(stats):
    res = "jac %.2f, recall %.2f, prec %.2f, rmsXY %.2f nm, rmsZ %.2f nm, avgX %.2f nm, avgY %.2f nm, avgZ %.2f nm, avgPhScale %.2f, avgPhDelta %.2f" % \
            (stats["jac"], stats["recall"], stats["precision"], 1e9 * stats["rmsXY"], \
             1e9 * stats["rmsZ"], 1e9 * stats["averageDx"], 1e9 * stats["averageDy"], \
             1e9 * stats["averageDz"], stats["avgPhotonsScale"], stats["avgPhotonsDelta"])   
    return res


#===============================================================================
# GroundTruthSettings
#===============================================================================

class GroundTruthSettings(ParamsBaseClass):
    """Class for matching ground-truth data.
            
    """
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.GroundTruths.GroundTruthHungarian")

        paramsThis = [
            Param("tolXY", 250e-9, 
                  friendlyName = "Tolerance XY", 
                  guiSiPrefix = True, 
                  guiSuffix = "m", 
                  guiStep = 10e-9),
                      
            Param("tolZ", 500e-9, 
                  friendlyName = "Tolerance Z", 
                  guiSiPrefix = True, 
                  guiSuffix = "m", 
                  guiStep = 10e-9)
            ]
        
        ParamsBaseClass.__init__(self, params + paramsThis, **kwargs)

    @property
    def name(self):
        return "Ground Truth"


class _GroundTruthStatsCommon():
    def SortedStatistics(self, locsToProcess, goodnessFunc, nPatch = 100, alreadySorted = False):
        self.logger.info("SortedStatistics...")
        
        # Initialize variables
        self.Clear()
        interpJacDict = {}
        jacs = []
        rmsXYs = []
        rmsZs = []
        maxDistXY = []
        
        if not alreadySorted:
            # Sort points
            locs = sorted(locsToProcess, key = goodnessFunc, reverse = True)
        else:
            locs = locsToProcess
            
        # Make patches
        nPatch = max(1, min(nPatch, len(locsToProcess)))
        locsPerPatch = float(len(locs)) / nPatch
        patchStarts = np.zeros((nPatch,), dtype = int)
        patchEnds = np.zeros((nPatch,), dtype = int)
        lastPatchEnd = 0
        for i in range(nPatch):
            # Compose slice for patch
            iS = lastPatchEnd
            iE = min(int(round((i + 1) * locsPerPatch)), len(locs))
            if i == nPatch - 1:
                iE = len(locs)
            patchStarts[i] = iS
            patchEnds[i] = iE
            lastPatchEnd = iE

        # Do statistics on ground-truth
        if self._fseries.HasGroundTruth():
            jacs = np.zeros((nPatch,))
            rmsXYs = np.zeros((nPatch,))
            rmsZs = np.zeros((nPatch,))
            precs = np.zeros((nPatch,))
            maxDistXY = np.zeros((nPatch,))
            
            for i in range(nPatch):
                # Compose slice for patch
                iS, iE = patchStarts[i], patchEnds[i]
                if iS >= iE:
                    # No element in patch
                    continue
                
                # Add locations from patch
                self.AddLocations(locs[iS:iE])
                s = self.GetStats()
                jacs[i] = s["jac"]
                rmsXYs[i] = s["rmsXY"]
                rmsZs[i] = s["rmsZ"] if s["rmsZ"] > 1e-16 else 1.0 
                precs[i] = s["precision"]
                maxDistXY[i] = max([l.distXY for l in locs[iS:iE]])
                
            # Save interpolated values
            extrapValue = 1.0
            jacInterp = np.concatenate((np.array([2.0, 4.0]), np.arange(5.0, 90.0 + 0.1, 5.0)))
            sl = slice(None, np.argmax(jacs) + 1)
            rmsXYInterp = np.interp(jacInterp, jacs[sl], rmsXYs[sl], left = extrapValue, right = extrapValue)
            rmsZInterp = np.interp(jacInterp, jacs[sl], rmsZs[sl], left = extrapValue, right = extrapValue)
            precInterp = np.interp(jacInterp, jacs[sl], precs[sl], left = extrapValue, right = extrapValue)
            
            for jac, trmsXY, trmsZ, tprec in zip(jacInterp, rmsXYInterp, rmsZInterp, precInterp):
                interpJacDict["rmsXY_jac_%.1f" % (jac)] = trmsXY
                interpJacDict["rmsZ_jac_%.1f" % (jac)] = trmsZ
                interpJacDict["prec_%.1f" % (jac)] = tprec
                
        sortingRes = {"jac": jacs, "rmsXY": rmsXYs, "rmsZ": rmsZs,
            "interpolated": interpJacDict, "maxDistXY": maxDistXY}
        self.logger.info("SortedStatistics done.")
        return locs, sortingRes
    
    def GetStatsStr(self):
        return GetStatsStr(self.GetStats())

#===============================================================================
# GroundTruthStatsPy
#===============================================================================

class GroundTruthStatsPy(object, _GroundTruthStatsCommon):

    def __init__(self, fseries, groundTruth, nrGroundTruthMolecules = None):
        self.logger = logging.getLogger("SMolPhot.GroundTruthStatsPy")
        self._fseries = fseries
        self._groundTruth = groundTruth
        self._nrGroundTruthMolecules = nrGroundTruthMolecules
        if nrGroundTruthMolecules is None:
            self._nrGroundTruthMolecules = self._fseries.GetNumberOfActMolecules()
        self.Clear()
        
    def MatchFrame(self, frame, locs):
        coordsArrayTest = GetLocsCoordsArray(locs)
        # CoordsTrue array    
        actIndices = frame._fseries._actFrameMap[frame.nr]
        coordsTrue = frame._fseries._actCoordsArray[actIndices, :]
        
        # Dist and cost matrices
        distsXY = spatial.distance.cdist(coordsArrayTest[:, :2], coordsTrue[:, :2])
        distsZ = spatial.distance.cdist(coordsArrayTest[:, 2:], coordsTrue[:, 2:])
        distsXY = np.where(distsXY > self._groundTruth.tolXY, 1.0, distsXY)
        distsZ = np.where(distsZ > self._groundTruth.tolZ, 1.0, distsZ)
        costMatrix = np.sqrt(distsXY ** 2.0 + distsZ ** 2.0)
        
        # Matching
        idsTest, _idsTrue = optimize.linear_sum_assignment(costMatrix)
        cost = costMatrix[idsTest, _idsTrue]
      
        # Filter out unsucessful matches
        resIndices = np.argwhere(cost < 1.0).ravel()
        if len(resIndices) <= 0:
            return np.array([]), np.array([]), np.array([]), np.array([])
        else:
            resIdsTest = idsTest[resIndices] 
            resIdsTrue = actIndices[_idsTrue[resIndices]]
            resDistsXY, resDistsZ = distsXY[resIdsTest, _idsTrue[resIndices]], distsZ[resIdsTest, _idsTrue[resIndices]]
            return resIdsTest, resIdsTrue, resDistsXY, resDistsZ
        
    def AddLocations(self, locs):
        self.logger.info("Add locations...")
        
        framesToMatch = set([loc.frameNr for loc in locs])
        # Remove stats from already match frames
        for frameNr in framesToMatch:
            if frameNr in self._matchingResults:
                # Frame already added previously, need to remove, match frame
                # again and then read the statistics
                frameLocs = self._locsByFrames[frameNr]
                frame = self._fseries.GetOriginalFrame(frameNr)
                self._RemoveFrameMatchingResults(frameLocs, frame, self._matchingResults[frameNr])
        
        # Update _locsByFrames dict
        for loc in locs:
            self._locsByFrames[loc.frameNr].append(loc)
        
        # Add stats for every frame
        for frameNr in framesToMatch:
            frameLocs = self._locsByFrames[frameNr]
            frame = self._fseries.GetOriginalFrame(frameNr)
            matchRes = self.MatchFrame(frame, frameLocs) 
            self._AddFrameMatchingResults(frameLocs, frame, matchRes)
            self._matchingResults[frameNr] = matchRes
            
        self.logger.info("Add locations done")
            
    def RemoveLocations(self, locsInitial):
        self.logger.info("Remove locations...")
        
        # Remove locations that is not in stats
        locs = []
        for loc in locsInitial:
            if not loc in self._locsByFrames[loc.frameNr]:
                continue
            locs.append(loc)
        
        # Remove stats from already match frames
        framesToMatch = set([loc.frameNr for loc in locs])
        for frameNr in framesToMatch:
            if frameNr in self._matchingResults:
                # Frame already added previously, need to remove, match frame
                # again and then read the statistics
                frameLocs = self._locsByFrames[frameNr]
                frame = self._fseries.GetOriginalFrame(frameNr)
                self._RemoveFrameMatchingResults(frameLocs, frame, self._matchingResults[frameNr])
        
        # Update _locsByFrames dict
        for loc in locs:
            self._locsByFrames[loc.frameNr].remove(loc)
            
        # Add stats for every frame
        for frameNr in framesToMatch:
            frameLocs = self._locsByFrames[frameNr]
            frame = self._fseries.GetOriginalFrame(frameNr)
            matchRes = self.MatchFrame(frame, frameLocs) 
            self._AddFrameMatchingResults(frameLocs, frame, matchRes)
            self._matchingResults[frameNr] = matchRes
            
        self.logger.info("Remove locations done")
            
    def Clear(self):
        self.logger.info("Clear...")
        self._locsByFrames = defaultdict(list)
        self._matchingResults = {}
        
        self.nTruePos = 0
        self.nFalsePos = 0
        self.nFalseNeg = 0
        self.distXYSqrSum = 0.0
        self.distZSqrSum = 0.0
        self.totalDx = 0.0
        self.totalDy = 0.0
        self.totalDz = 0.0
        self.photonsScale = 0.0
        self.photonsDelta = 0.0
        self.logger.info("Clear done")

            
    def GetStats(self):
        self.logger.info("Get stats...")
        nTP, nFP, nFN = self.nTruePos, self.nFalsePos, self.nFalseNeg
        
        if nTP > 0:
            # Need to add actMolsNotProcessed to nFN
            actMolsProcessed = nFN + nTP
            actMolsNotProcessed = self._nrGroundTruthMolecules - actMolsProcessed
            nFN += actMolsNotProcessed
            
            rmsXY = float(np.sqrt(self.distXYSqrSum / nTP)) 
            rmsZ = float(np.sqrt(self.distZSqrSum / nTP)) 
            jac = 100.0 * float(nTP) / (nFN + nFP + nTP)
            recall = float(nTP) / (nFN + nTP)
            precision = float(nTP) / (nFP + nTP)
            averageDx = self.totalDx / nTP
            averageDy = self.totalDy / nTP
            averageDz = self.totalDz / nTP
            avgPhotonsScale = self.photonsScale / nTP
            avgPhotonsDelta = self.photonsDelta / nTP
            
        else:
            rmsXY = None
            rmsZ = None
            jac = None
            recall = None
            precision = None
            averageDx = None 
            averageDy = None
            averageDz = None
            avgPhotonsScale = None
            avgPhotonsDelta = None
            
            
        res = {"rmsXY": rmsXY,
                 "rmsZ": rmsZ,
                 "jac": jac,
                 "recall": recall,
                 "precision": precision,
                 "nTP": nTP,
                 "nFP": nFP,
                 "nFN": nFN,
                 "averageDx": averageDx,
                 "averageDy": averageDy,
                 "averageDz": averageDz,
                 "avgPhotonsScale": avgPhotonsScale,
                 "avgPhotonsDelta": avgPhotonsDelta
                 }
        self.logger.info("Get stats done")
        return res
    
    def _AddFrameMatchingResults(self, frameLocs, frame, matchRes, coef = 1.0):
        # coef = 1: add, coef = -1: remove
        idsTest, idsTrue, distsXY, distsZ = matchRes
            
        # Calculate average XY deviation
        for idTest, idTrue in zip(idsTest, idsTrue):
            testX, testY, testZ = frameLocs[idTest].coord
            trueX = self._fseries._actCoords[0][idTrue]
            trueY = self._fseries._actCoords[1][idTrue]
            trueZ = self._fseries._actCoords[2][idTrue]
            
            self.totalDx += coef * float(testX - trueX)
            self.totalDy += coef * float(testY - trueY)
            self.totalDz += coef * float(testZ - trueZ) if testZ is not None else 0.0
            # print "framenr: " , frameNr , " dx: " , float(testX - trueX) , " dy: " ,float(testY - trueY)
            
            # reported indensity deviation
            photonsTest = frameLocs[idTest].photons
            photonsTrue = self._fseries._actIs[idTrue]
            self.photonsScale += coef * float(photonsTest / photonsTrue)
            self.photonsDelta += coef * float(photonsTest - photonsTrue)
            
  
        # Optimal goodness
        # Remove
        for loc in frameLocs:
            loc.distXY = self._groundTruth.tolXY
            loc.distZ = self._groundTruth.tolZ
        # Add
        if coef > 0.0:
            for idTest, distXY, distZ in zip(idsTest, distsXY, distsZ): 
                frameLocs[idTest].distXY = distXY
                frameLocs[idTest].distZ = distZ
                

        nTP = len(distsXY)
        nFP = len(frameLocs) - nTP
        nFN = len(self._fseries._actFrameMap[frame.nr]) - nTP
        
        self.nTruePos += int(coef) * nTP
        self.nFalsePos += int(coef) * nFP
        self.nFalseNeg += int(coef) * nFN
        self.distXYSqrSum += coef * np.sum(distsXY ** 2.0)
        self.distZSqrSum += coef * np.sum(distsZ ** 2.0)
        
    def _RemoveFrameMatchingResults(self, *args, **kwargs):
        self.logger.info("_RemoveFrameMatchingResults")
        kwargs.pop("coef", None)
        self._AddFrameMatchingResults(*args, coef = -1.0, **kwargs)


#===============================================================================
# GroundTruthStatsLibrary
#===============================================================================

try:
    from GroundTruthStatsCpp import GroundTruthStatsCpp
    from SMolPhot.PSFs import LocsToArrays
    
    class GroundTruthStatsLibrary(GroundTruthStatsCpp, _GroundTruthStatsCommon):
        
        def __init__(self, fseries, groundTruth, nrGroundTruthMolecules = None):
            self.logger = logging.getLogger("SMolPhot.GroundTruthStatsLibrary")
            self._fseries = fseries
            self._groundTruth = groundTruth
            if nrGroundTruthMolecules is None:
                nrGroundTruthMolecules = self._fseries.GetNumberOfActMolecules()
            self._nrGroundTruthMolecules = nrGroundTruthMolecules
        
            GroundTruthStatsCpp.__init__(self, self._fseries._actFrameNr, \
                                     self._fseries._actCoords, \
                                     self._fseries._actIs, \
                                     nrGroundTruthMolecules, \
                                     self._groundTruth.tolXY, \
                                     self._groundTruth.tolZ)
            
        def MatchFrame(self, frame, locs):
            coordsArrayTest = GetLocsCoordsArray(locs)
            frameNrs = np.ascontiguousarray(np.squeeze([loc.frameNr for loc in locs]), dtype = int)
            coords = (np.ascontiguousarray(coordsArrayTest[:, 0]), 
                      np.ascontiguousarray(coordsArrayTest[:, 1]),
                      np.ascontiguousarray(coordsArrayTest[:, 2]))
            photons = np.ascontiguousarray(np.squeeze([loc.photons for loc in locs]), dtype = float)
            return GroundTruthStatsCpp.MatchFrame(self, frame.nr, frameNrs, coords, photons)
    
        def AddLocations(self, locs):
            self.logger.info("Add locations (%d)" % (len(locs)))
            distsXY, distsZ = GroundTruthStatsCpp.AddLocations(self, *LocsToArrays(locs))
            for i, loc in enumerate(locs):
                loc.distXY = distsXY[i]
                loc.distZ = distsZ[i]
            
        def RemoveLocations(self, locs):
            self.logger.info("Remove locations (%d)" % (len(locs)))
            GroundTruthStatsCpp.RemoveLocations(self, *LocsToArrays(locs))
            
        def Clear(self):
            self.logger.info("Clear...")
            GroundTruthStatsCpp.Clear(self)
            
        def GetStats(self):
            self.logger.info("Get stats...")
            return GroundTruthStatsCpp.GetStats(self)

    GroundTruthStats = GroundTruthStatsLibrary
except ImportError:
    print "Unable to import c++ ground truth, using pure python"
    GroundTruthStats = GroundTruthStatsPy
    
if __name__ == "__main__":
    pass
