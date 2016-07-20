"""Implements simple local maxima localizer.

"""
import logging
import numpy as np
from collections import OrderedDict
from SMolPhot.Components.BaseClasses import Param
from _Common import LocalizerBase, AddLocsWoDuplicates, CalcMinimumDistance

__all__ = ["LocalMaximaLocalizer"]

class LocalMaximaLocalizer(LocalizerBase):
    """Simple local maxima localizer.
    
    Params:
        detThreshold (float): threshold = mean + detThreshold * std
        fitPixels (int): half-width of the rectangular fit region
    
    """
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Localizers.LocalMaximaLocalizer")
    
        self._initialFitLocs = OrderedDict([("Local maxima center", "localMaxima"),
                               ("Weighted centroid", "weightedCentroid")])
        
        paramsThis = [Param("detThreshold", \
                            2000.0, \
                            friendlyName = "Maxima threshold", \
                            guiDecimals = 4, \
                            guiSuffix = " ph", \
                            guiStep = 1.0, \
                            guiLimits = (1.0, 100000)),
                      
                      Param("fitPixels", \
                            1, \
                            friendlyName = "Region half-width", \
                            guiSuffix = " px", \
                            guiLimits = (1, 10)),
                      
                      Param("_minFitPixels", \
                            1, \
                            friendlyName = "Min region half-width", \
                            guiSuffix = " px", \
                            guiLimits = (1, 10)),
      
                      Param("_initialMaxDist", \
                            200e-9, \
                            friendlyName = "Initial max dist", \
                            guiSiPrefix = True, \
                            guiSuffix = "m", \
                            guiStep = 10e-9),
                      
                      Param("_borderRegionWidth", \
                            300e-9, \
                            friendlyName = "Border region", \
                            guiSiPrefix = True, \
                            guiSuffix = "m", \
                            guiStep = 10e-9,
                            guiLimits = (0.0, 1.0)),
                      
                        Param("_duplicateMinDist", \
                            200e-9, \
                            friendlyName = "Duplicate min dist", \
                            guiSiPrefix = True, \
                            guiSuffix = "m", \
                            guiStep = 10.0e-9,
                            guiLimits = (0.0, 1000e-9)),

                        Param("_minDistSimFit", \
                            1000e-9, \
                            friendlyName = "Simultaneous fit min dist", \
                            guiSiPrefix = True, \
                            guiSuffix = "m", \
                            guiStep = 10.0e-9,
                            guiLimits = (0.0, 10000e-9)),
                      
                      
                        Param("_minGoodness", \
                            0.5, \
                            friendlyName = "Min goodness", \
                            guiDecimals = 5, \
                            guiStep = 0.01,
                            guiLimits = (0.0, 1e12)),

                      Param("_initialFitLoc",
                            "localMaxima",
                            friendlyName = "Initial fit location",
                            guiType = "list",
                            guiValues = self._initialFitLocs
                            ),
                      
                      Param("_minOffset", \
                            -20.0, \
                            friendlyName = "Min offset", \
                            guiDecimals = 4, \
                            guiSuffix = " ph", \
                            guiStep = 1.0, \
                            guiLimits = (-1000, 1000)),
                      
                    Param("_noiselevel", \
                            10.0, \
                            friendlyName = "Noise level", \
                            guiDecimals = 4, \
                            guiSuffix = " ph", \
                            guiStep = 1.0, \
                            guiLimits = (-1000, 1000)),
                      
                    Param("_minArea", \
                            4, \
                            friendlyName = "Min area", \
                            guiSuffix = " px", \
                            guiLimits = (1, 10)),
                      ]
        
        LocalizerBase.__init__(self, params + paramsThis, **kwargs)
    
    def OptimizeForSeries(self, fseries):
        pass
    
    def FindMolecules(self, inputFrame, psf, axialCalibrator, fitMode = "z"):
        
        def _FitMolecule(frame, (x, y), fitPixels):
            # Fit            
            locFit = psf.Fit(frame, (x, y), fitPixels, \
                             fitMode = fitMode, \
                             initialMaxDist = self._initialMaxDist)
            
            if locFit is None:
                return None
                        
            if min(locFit.x, locFit.y) < self._borderRegionWidth or \
                locFit.x + self._borderRegionWidth > inputFrame.sizeX or \
                locFit.y + self._borderRegionWidth > inputFrame.sizeY:
                self.logger.debug("Dismissed, in border requion")
                return None
            
            if locFit.z is not None:
                if locFit.z < axialCalibrator._calibFrom or locFit.z > axialCalibrator._calibTo:
                    self.logger.debug("Dismissed, outside axial calibration range %.2f nm" % (1e9 * locFit.z))
                    return None

            if locFit.goodnessValue < self._minGoodness:
                self.logger.debug("Dismissed, to bad goodness %.2f" % (locFit.goodnessValue))
                return None
            
            if locFit.offset < self._minOffset:
                self.logger.debug("Dismissed, to bad offset %.2f" % (locFit.offset))
                return None
            
            return locFit
        
        def FindLocalMaximas(frame, threshold, mode, noiseLevel = -np.inf, minArea = 1):
            maxIndicesInitial = np.where(frame.data > threshold)
            pixelsX, pixelsY = frame.pixelsX, frame.pixelsY
            
            maxIndicesX, maxIndicesY = [], []
            for xI, yI in zip(maxIndicesInitial[0], maxIndicesInitial[1]):
                #frame.logger.debug("Potential local maxuma at pixel (%d, %d)" % (xI, yI))
                if min(xI, yI) < 1 or \
                    xI + 1 >= pixelsX or \
                    yI + 1 >= pixelsY:
                    #frame.logger.debug("Dismissed, pixel on edge")
                    continue
                
                pixelValue = frame.data[xI, yI]
                directlyAroundData = frame.data[(xI - 1):(xI + 2), (yI - 1):(yI + 2)]
                maxAround = np.max(directlyAroundData)
                #minAround = np.min(directlyAroundData)
                area = (directlyAroundData > noiseLevel).sum()
                
                if maxAround > pixelValue:
                    # Not local maxima
                    #frame.logger.debug("Dismissed, not local maxima")
                    continue
                
                if area < minArea:
                    # area to small
                    continue
                
                if mode == "localMaxima":
                    maxIndicesX.append(xI)
                    maxIndicesY.append(yI)
                elif mode == "weightedCentroid":
                    tmpX, tmpY = np.meshgrid(np.arange(directlyAroundData.shape[0]), \
                                             np.arange(directlyAroundData.shape[1]), \
                                             indexing = "ij")
                    centroidXI = np.average(tmpX, weights = directlyAroundData) - 1.0 + xI
                    centroidYI = np.average(tmpY, weights = directlyAroundData) - 1.0 + yI
                    maxIndicesX.append(centroidXI)
                    maxIndicesY.append(centroidYI)
                    #print "centroid vs local maximum", (xI, yI), (centroidXI, centroidYI)
                else:
                    raise NotImplementedError()
                    
                    
            maxIndices = (np.array(maxIndicesX), np.array(maxIndicesY))
            maxCoords = frame.GetPixelCoords(maxIndices)
            return maxIndices, maxCoords
        
        class FakeRegion(object):
            pass
        
        def RecursiveSearch(iteration, frame, fitPixels, locsAlreadyFound = []):
            if fitPixels < self._minFitPixels:
                return locsAlreadyFound
            
            iterationInfo = {"iteration": iteration, "fitPixels": fitPixels, \
                             "startFrame": frame}
            self._iterationsInfoTmp.append(iterationInfo)
            
            self.logger.info("RecursiveSearch iteration %d, fitPixels %d" % (iteration, fitPixels))
            maxIndices, maxCoords = FindLocalMaximas(frame, self.detThreshold, \
                                                           mode = self._initialFitLoc, \
                                                           noiseLevel = self._noiselevel,
                                                           minArea = self._minArea)
            
            if len(maxIndices[0]) == 0:
                # No local maximums
                return locsAlreadyFound
            
            regions = []
            locsThis = []
            for x, y, xI, yI in zip(maxCoords[0], maxCoords[1], maxIndices[0], maxIndices[1]):                
                fakeRegion = FakeRegion()
                fakeRegion.bbox = (-1, -1, 2, 2)
                fakeRegion.offset = np.array([int(round(xI)), int(round(yI))])
                regions.append(fakeRegion)
                
                locFit = _FitMolecule(frame, (x, y), fitPixels)
                
                if locFit is None:
                    continue                    
                
                locFit.iteration = iteration
                locsThis.append(locFit)
                
            locsThis.sort(key = lambda loc: -loc.goodnessValue)
            locs, locsThis = AddLocsWoDuplicates(locsAlreadyFound, locsThis, \
                                                 self._minDistSimFit, \
                                                 self._duplicateMinDist, \
                                                 self.logger)                
            
            # Remove found PSFs
            frameNext = frame.CopyFrame().SubtractLocs(psf, locsThis)
            
            # Save additional iterationInfo
            iterationInfo["psfsToDelete"] = frame.data - frameNext.data
            iterationInfo["frameNext"] = frameNext
            iterationInfo["locs"] = locsThis
            iterationInfo["regions"] = []
            iterationInfo["internalRegions"] = regions
            
            self.logger.info("Locs from recursiveSearch iteration %d: %d" % (iteration, len(locsThis)))
            
            fitPixelsNext = fitPixels - 1 if len(locsThis) == 0 else fitPixels
            locs = RecursiveSearch(iteration + 1, frameNext, fitPixelsNext, locsAlreadyFound = locs)
            
            return locs
        
        # Main ----------------------------------------------------------------- 

        self._iterationsInfoTmp = []
        self.logger.info("Find molecules frame #%d, fitmode %s" % (inputFrame.nr, fitMode))
        locs = RecursiveSearch(0, inputFrame, self.fitPixels)
        locs = CalcMinimumDistance(locs)
        return locs
        
    @property
    def name(self):
        return "Local Maxima Localizer"

if __name__ == '__main__':
    pass
