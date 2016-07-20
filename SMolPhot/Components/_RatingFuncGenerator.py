import os
import cPickle
import itertools
import numpy as np
# import cProfile, pstats

from os import path
from time import time
from collections import defaultdict
from SMolPhot.Components import GroundTruthStats
from SMolPhot.PSFs import LocsToArrays


__all__ = []

class RatingFuncGenerator(object):
    def __init__(self, fseries, groundTruth, locs):
        self._fseries = fseries
        self._groundTruth = groundTruth
        self._locs = locs
        
        self._nrToCombine = 3
        self._powersToTest = np.arange(-3.0, 3.0 + 1e-9, 1.0)
        self._listOfProperties = ["lsError", "photons", "amp", "absz", "sigmaX", \
                                 "sigmaY", "absoffset", "fitStdAmp", "fitStdX0", \
                                 "fitStdY0", "fitStdZ0", "fitStdOffset", \
                                 "minFitDistXY", "nrOfPixels", \
                                 "distFitWeightedCentroid", "area"]
        
        #self._listOfProperties = ["lsError", "photons", "nrOfPixels", "fitStdAmp", \
        #                          "fitStdOffset"]
        
    def GetLocsProperties(self, locs):
        res = {}
        for propertyName in self._listOfProperties:
            array = np.squeeze(np.array([getattr(loc, propertyName) for loc in locs], dtype = float))
            res[propertyName] = array
        return res
    
    def GetErrorFuncStats(self, errorFuncData):
        # Import C++ GroundTruth
        from GroundTruthStatsCpp import GroundTruthStatsCpp
        
        # Add every point to GroundTruthStats
        gtStat = GroundTruthStatsCpp(self._fseries._actFrameNr, \
                                     self._fseries._actCoords, \
                                     self._fseries._actIs, \
                                     self._fseries.GetNumberOfActMolecules(), \
                                     self._groundTruth.tolXY, \
                                     self._groundTruth.tolZ)
        gtStat.AddLocations(*LocsToArrays(self._locs))
        
        # Save initial
        res = defaultdict(list)
        res["partRemoved"].append(0)
        for k, v in gtStat.GetStats().iteritems():
            res[k].append(v)
    
        # Start removing points based on errorFuncData
        locs = self._locs
        nrOfAlreadyRemoved = 0
        for part, errorFunc in errorFuncData:
            # Sort point based on errorFunc
            locs.sort(key = errorFunc, reverse = True)
            
            # Remove points
            nrOfLocsToRemove = int(round(part * len(self._locs))) - nrOfAlreadyRemoved 
            locsToRemove = locs[:nrOfLocsToRemove]
            gtStat.RemoveLocations(*LocsToArrays(locsToRemove))
            stats = gtStat.GetStats()
            
            # Save
            res["partRemoved"].append(part)
            for k, v in stats.iteritems():
                res[k].append(v)
            
            # newLocs
            toRemoveSet = set(locsToRemove)
            locsNew = [loc for loc in locs if loc not in toRemoveSet]
            
            # Prepare for next loop
            nrOfAlreadyRemoved += nrOfLocsToRemove
            locs = locsNew
        
        # Convert res to numpy arrays
        for k, v in res.iteritems():
            res[k] = np.array(v)
        return res
    
    def GenerateGoodnessFunctions(self, saveFolder):
        res = []
        locs = self._locs
        for comb in itertools.combinations(self._listOfProperties, self._nrToCombine):
            for pwrs in itertools.product(*([self._powersToTest.tolist()] * self._nrToCombine)):
                goodnessFuncTmp = lambda loc: np.prod([getattr(loc, prop) ** pwr for prop, pwr in zip(comb, pwrs)])
                
                statTmp = GroundTruthStats(self._fseries, self._groundTruth)
                _, sortingResTmp = statTmp.SortedStatistics(locs, goodnessFuncTmp)
                rmsXY_jac_20 = sortingResTmp["interpolated"]["rmsXY_jac_20.0"]
                print comb, pwrs, "rmsXY_jac_20", 1e9 * rmsXY_jac_20
                
                res.append((comb, pwrs, sortingResTmp))
       
        res.sort(key = lambda x: np.sum(x[2]["rmsXY"]))
       
        # Save to file
        if not path.isdir(saveFolder):
            os.mkdir(saveFolder)
        saveFile = path.join(saveFolder, "goodnesTest_%d.csv" % (self._nrToCombine))
        with open(saveFile, "w") as f:
            first = True
            for comb, pwrs, sortingRes in res[:100000]:
                if first:
                    sortedReskeys = sorted(sortingRes["interpolated"].keys())
                    f.write("Comb;Pwrs;")
                    f.write(";".join(sortedReskeys))
                    f.write("\n")
                    first = False
                    
                f.write("%s;%s;" % (comb, pwrs))
                for k in sortedReskeys: 
                    interpV = sortingRes["interpolated"][k]
                    coef = 1e9
                    f.write("%.4f;" % (coef * interpV))
                f.write("\n")

    def GenerateOptimalErrorFunction(self, partToRemove, nrOfIterationsList, bestnessCriteriaFunc, saveFolder):
        startTime = time()
        allInfo = {}
        txtInfo = []

        for nrOfIterations in nrOfIterationsList:
            iterationInfo = self._GenerateErrorFunction(partToRemove, nrOfIterations, bestnessCriteriaFunc)
            
            # Info
            allInfo[nrOfIterations] = iterationInfo
            txtInfo.append("nrOfIterations %d" % (nrOfIterations))
            for part, bestComb, bestPwrs, bestStat in zip(iterationInfo["partRemoved"], iterationInfo["bestComb"], iterationInfo["bestPwrs"], iterationInfo["bestStats"]):
                txtInfo.append("%f %s %s, jac %.4f %.4f nm, %.4f nm" % \
                    (part, str(bestComb), str(bestPwrs), bestStat["jac"], \
                     1e9 * bestStat["rmsXY"], 1e9 * bestStat["rmsZ"]))
                
            txtInfo.append("python code:")
            txtInfo.append('errorFunctions["%s_%d"] = [' % (self._fseries.name, nrOfIterations))
            for part, bestComb, bestPwrs in zip(iterationInfo["partRemoved"], iterationInfo["bestComb"], iterationInfo["bestPwrs"]):
                if part <= 0.0:
                    continue
                funcText = " * ".join(["loc.%s ** %.1f" % (param, pwr) for param, pwr in zip(bestComb, bestPwrs)])
                txtInfo.append("(%.6f, lambda loc: %s)," % (part, funcText))
            txtInfo.append("]")
             
            txtInfo.append("")
            
        # Save all data
        if not path.isdir(saveFolder):
            os.mkdir(saveFolder)
        with open(path.join(saveFolder, "data.pkl"), "wb") as output:
            cPickle.dump(iterationInfo, output, protocol = -1)
       
        # Save txt
        with open(path.join(saveFolder, "info.txt"), "w") as output:
            output.write("\n".join(txtInfo))
        
        totalTime = time() - startTime
        print "Total time %.2f h" % (totalTime / 3600.0)
        return allInfo
       
    def _GenerateErrorFunction(self, totalPartToRemove, nrOfIterations, bestnessCriteriaFunc):
        
        def AddToIterationInfo(partToRemove, iteration, comb, pwrs, stats, alternativeFuncs = None):
            iterationInfo["iteration"].append(iteration)
            iterationInfo["partRemoved"].append(partToRemove)
            iterationInfo["bestComb"].append(comb)
            iterationInfo["bestPwrs"].append(pwrs)
            iterationInfo["bestStats"].append(stats)
            iterationInfo["alternativeFuncs"].append(alternativeFuncs)
            for k, v in stats.iteritems():
                iterationInfo[k].append(v)
        
        startTime = time()
        iterationInfo = defaultdict(list)
        
        locs = self._locs
        totLocsToRemove = int(round(totalPartToRemove * len(self._locs)))
        nrOfAlreadyRemoved = 0
        print "GenerateErrorFunction: Remove %d locs in %d iterations" % (totLocsToRemove, nrOfIterations)
        for iteration in range(nrOfIterations):
            print "iteration %d (# of locs %d)" % (iteration, len(locs))
            # Find best error functions to remove nrOfLocsToRemove
            #nrOfLocsToRemove = int(round((iteration + 1) * float(totLocsToRemove))) / nrOfIterations - alreadyRemoved
            partToRemove = (iteration + 1) / float(nrOfIterations) * totalPartToRemove
            nrOfPointsToRemove = int(round(partToRemove * float(len(self._locs)))) - nrOfAlreadyRemoved
            (bestComb, bestPwrs, bestStats), newLocs, initialStats, alternativeFuncs = \
                self._FindBestFunctionIteration(locs, nrOfPointsToRemove, bestnessCriteriaFunc)
            nrOfAlreadyRemoved += nrOfPointsToRemove
            print
            
            # Save iteration info
            if iteration == 0:
                AddToIterationInfo(0.0, -1, bestComb, bestPwrs, initialStats)
            AddToIterationInfo(partToRemove, iteration, bestComb, bestPwrs, bestStats, alternativeFuncs)
            
            # Save newLocs for next iteration
            locs = newLocs
        
        # Convert to float
        toArrayList = ["partRemoved", "jac", "rmsXY", "rmsZ", "recall"]
        for k in toArrayList:
            iterationInfo[k] = np.array(iterationInfo[k])
        
        totalTime = time() - startTime
        print "Total iteration time %.2f h" % (totalTime / 3600.0)
        
        return iterationInfo
    
    def _FindBestFunctionIteration(self, locs, nrOfLocsToRemove, bestnessCriteriaFunc):
        startTime = time()
        nrOfLocsToRemove = min(nrOfLocsToRemove, len(locs))
        
        # Import C++ GroundTruth
        from GroundTruthStatsCpp import GroundTruthStatsCpp
        # profile = cProfile.Profile()
        # profile.enable()
        
        # Generate list of properties
        properties = self.GetLocsProperties(locs)
        
        # Add every point to GroundTruthStats
        gtStat = GroundTruthStatsCpp(self._fseries._actFrameNr, \
                                     self._fseries._actCoords, \
                                     self._fseries._actIs, \
                                     self._fseries.GetNumberOfActMolecules(), \
                                     self._groundTruth.tolXY, \
                                     self._groundTruth.tolZ)
        gtStat.AddLocations(*LocsToArrays(locs))
        
        # Generate combinations
        res = []
        nrCombinations = 0
        for comb in itertools.combinations(self._listOfProperties, self._nrToCombine):
            for pwrs in itertools.product(*([self._powersToTest.tolist()] * self._nrToCombine)):
                # print comb, pwrs
                # Calc error values and sort
                errorValues = np.prod([properties[prop] ** pwr for prop, pwr in zip(comb, pwrs)], axis = 0)                
                sortedIndices = np.argsort(-errorValues)
                locsToRemove = [locs[i] for i in sortedIndices[:nrOfLocsToRemove]]

                # Remove nrOfLocsToRemove most bad points
                locsArrays = LocsToArrays(locsToRemove)
                gtStat.RemoveLocations(*locsArrays)
                
                # Calculate statisticsa and add points back
                stats = gtStat.GetStats()
                gtStat.AddLocations(*locsArrays)
                
                # Save result
                res.append((comb, pwrs, stats, locsToRemove))
                nrCombinations += 1
        
        # Sort all combinations by largest jac and by smallest rmsXY
        res.sort(key = bestnessCriteriaFunc)  # lambda x: (-x[2]["jac"], x[2]["rmsXY"])
        # res.sort(key = lambda x: (x[2]["rmsXY"]))
        bestComb, bestPwrs, bestStats, bestLocsToRemove = res[0]
        alternativeFuncs = res[:50]
        
        # Generate new locs array
        toRemove = set(bestLocsToRemove)
        locsNew = [loc for loc in locs if loc not in toRemove]
        
        # Print info
        totalTime = time() - startTime
        estimateTimePerM = totalTime / nrCombinations * 1e6
        initialStats = gtStat.GetStats()
        print "initial jac %.4f, rmsXY %.4f, rmsZ %.4f" % (initialStats["jac"], \
                                                           1e9 * initialStats["rmsXY"], \
                                                           1e9 * initialStats["rmsZ"]) 
        
        print "best jac %.4f, rmsXY %.4f, rmsZ %.4f" % (bestStats["jac"], \
                                                           1e9 * bestStats["rmsXY"], \
                                                           1e9 * bestStats["rmsZ"]) 
        print "best func", bestComb, bestPwrs
        print "nr of combinations %d, time elapsed %.1f s, estimeted time per 1M %.2f h" % \
            (nrCombinations, totalTime, estimateTimePerM / 3600.0)
 
        return (bestComb, bestPwrs, bestStats), locsNew, initialStats, alternativeFuncs
        
if __name__ == "__main__":
    pass
