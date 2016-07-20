import logging
import numpy as np
# import cProfile, pstats

from collections import OrderedDict
from SMolPhot.Components.BaseClasses import Param
from _Common import LocalizerBase, CalcMinimumDistance, SimpleRegion

__all__ = ["IterativeLocalizer"]

#===============================================================================
# IterativeLocalizer
#===============================================================================

class IterativeLocalizer(LocalizerBase):
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Localizers.IterativeLocalizer")
    
        self._potentialLocsMode = OrderedDict([("Local maxima", "localMaxima"),
                                               ("Above treshold", "aboveTreshold")])
        
        self._initialFitLocs = OrderedDict([("Local maxima center", "localMaxima"),
                               ("Weighted centroid", "weightedCentroid")])
        
        self._fitPixelModes = OrderedDict([("Reduce min max", "reduceMinMax"),
                               ("Reduce min", "reduceMin"),
                               ("All", "all")])
        
        paramsThis = [
            # General params
                      
            Param("_saveIterationInfo", False,
                  friendlyName = "Iteration info"),
            
            Param("_maxIterations", 10,
                  friendlyName = "Max iterations",
                  guiLimits = (1, 100)),
                      
            # Find potential locations params
            
            Param("_potentialLocMode", "localMaxima",
                  friendlyName = "Potential loc mode",
                  guiType = "list",
                  guiValues = self._potentialLocsMode),
              
            Param("_detThreshold", 50.0,
                  friendlyName = "Detection threshold",
                  guiDecimals = 4,
                  guiSuffix = " ph",
                  guiStep = 1.0,
                  guiLimits = (1.0, 100000)),
              
            # Potential locations filtering params
            
            Param("_noiselevel", 40.0,
                  friendlyName = "Noise level",
                  guiDecimals = 4,
                  guiSuffix = " ph",
                  guiStep = 1.0,
                  guiLimits = (-1000, 1000)),
              
            Param("_minArea", 4,
                  friendlyName = "Min area",
                  guiSuffix = " px",
                  guiLimits = (1, 10)),
              
            # Fitting params

            Param("_minFitPixels", 3,
                  friendlyName = "Min region half-width",
                  guiSuffix = " px",
                  guiLimits = (1, 10)),
                      
            Param("_maxFitPixels", 5,
                  friendlyName = "Max region half-width",
                  guiSuffix = " px",
                  guiLimits = (1, 10)),
                      
            Param("_fitPixelMode", "reduceMinMax",
                  friendlyName = "Region half-width mode",
                  guiType = "list",
                  guiValues = self._fitPixelModes),
            
            Param("_initialMaxDist", 250e-9,
                  friendlyName = "Initial max dist",
                  guiSiPrefix = True,
                  guiSuffix = "m",
                  guiStep = 10e-9),
              
            Param("_borderRegionWidth", 300e-9,
                  friendlyName = "Border region",
                  guiSiPrefix = True,
                  guiSuffix = "m",
                  guiStep = 10e-9,
                  guiLimits = (0.0, 1.0)),
                      
            Param("_initialFitLoc", "localMaxima",
                  friendlyName = "Initial fit location",
                  guiType = "list",
                  guiValues = self._initialFitLocs),
              
            Param("_minOffset", -100.0,
                  friendlyName = "Min offset",
                  guiDecimals = 4,
                  guiSuffix = " ph",
                  guiStep = 1.0,
                  guiLimits = (-1000, 1000)),
              
            Param("_duplicateMinDist", 400e-9,
                  friendlyName = "Duplicate min dist",
                  guiSiPrefix = True,
                  guiSuffix = "m",
                  guiStep = 10.0e-9,
                  guiLimits = (0.0, 1000e-9)),
              
            Param("_minGoodness", 0.5,
                  friendlyName = "Min goodness",
                  guiDecimals = 5,
                  guiStep = 0.01,
                  guiLimits = (0.0, 1e99)),

        ]
        
        LocalizerBase.__init__(self, params + paramsThis, **kwargs)
    
    def OptimizeForSeries(self, fseries):
        pass

    def FindMolecules(self, inputFrame, psf, axialCalibrator, fitMode = "z"):
        self.logger.info("FindMolecules...")
        
        # Init variables
        frame = inputFrame.CopyFrame()
        self._iterationInfo = {}
        self._iterationsInfoTmp = []
        
        # Init potential locations collection
        potentialLocs = PotentialLocCollection(self, frame, psf, axialCalibrator, fitMode, \
                                               saveIterationInfo = self._saveIterationInfo)
        if self._fitPixelMode == "all":
            potentialLocs.SetFitPixelsRange(self._minFitPixels, self._maxFitPixels)
            
        # Do iterations
        for iteration in range(self._maxIterations):
            self.logger.info("FindMolecules iteration: %d" % (iteration))
            
            if self._saveIterationInfo:
                self._iterationInfo["iteration"] = iteration
                self._iterationInfo["startFrame"] = frame.CopyFrame()
                self._iterationsInfoTmp.append(self._iterationInfo)
                
            # Try to find on enough good locFit
            while True:
                locFit = potentialLocs.GetBestLocFit()
                if locFit is not None:
                    break
                
                # No location was found, need to loosen the parameters or to terminate
                if self._fitPixelMode == "reduceMinMax":
                    if potentialLocs.minFitPixels > self._minFitPixels:
                        potentialLocs.SetFitPixelsRange(potentialLocs.minFitPixels - 1, \
                                                        potentialLocs.minFitPixels - 1)
                    else:
                        break
                elif self._fitPixelMode == "reduceMin":
                    if potentialLocs.minFitPixels > self._minFitPixels:
                        potentialLocs.SetFitPixelsRange(potentialLocs.minFitPixels - 1, \
                                                        potentialLocs.maxFitPixels)
                    else:
                        break
                else:
                    break
            
            # No good locFit was found, terminate        
            if locFit is None:
                break
            
            # Remove found molecule and start over
            potentialLocs.AddLocationAndUpdate(locFit)
            self._iterationInfo = {}
                  
        # Warn if all iterations were used
        if iteration + 1 >= self._maxIterations:
            self.logger.warn("Maximum number of iterations (%d) exceeded at frame %d" % \
                             (self._maxIterations, frame.nr + 1))
                        
        # Calculate minFitDistXY for every loc
        locs = CalcMinimumDistance(potentialLocs.resultLocs)
        
        self.logger.info("FindMolecules done.")
        return locs
        
    @property
    def name(self):
        return "Iterative Localizer"


#===============================================================================
# PotentialLocCollection
#===============================================================================

class PotentialLocCollection(object):
    """This class holds all potential fitting locations (insede class PotentalLoc)
    and mainly provides GetBestFitLoc() and AddLocationAndUpdate() functionality.
    
    """
    def __init__(self, localizer, frame, psf, axialCalibrator, fitMode, saveIterationInfo = False):
        self.logger = localizer.logger
        self._localizer = localizer
        self._frame = frame
        self._psf = psf
        self._axialCalibrator = axialCalibrator
        self._fitMode = fitMode
        self._minFitPixels = localizer._maxFitPixels
        self._maxFitPixels = localizer._maxFitPixels
        self.saveIterationInfo = saveIterationInfo
        
        # Init
        self._potentialLocs = []
        self.resultLocs = []
        self.AddPotentialLocs()
        
    def AddPotentialLocs(self, bbox = None):
        self.logger.info("AddPotentialLocs: %s" % (str(bbox)))
        
        if self._localizer._potentialLocMode == "localMaxima":
            initialCoords = self._frame.FindMaximas(self._localizer._detThreshold, \
                                                    bbox = bbox, localMaxima = True)
        elif self._localizer._potentialLocMode == "aboveTreshold":
            initialCoords = self._frame.FindMaximas(self._localizer._detThreshold, \
                                                    bbox = bbox, localMaxima = False)
        else: 
            raise NotImplementedError()
        
        # Add new points that match conditions
        toAdd = []
        for coord in initialCoords:
            potentialLoc = PotentialLoc(self, coord)
            if potentialLoc.area < self._localizer._minArea:
                continue
            toAdd.append(potentialLoc)
        self._potentialLocs += toAdd
        
        # Add iteration info
        if self.saveIterationInfo:
            iterationInfo = self._localizer._iterationInfo
            iterationInfo["updateRegions"] = [SimpleRegion((0.0, 0.0), bbox)]
            iterationInfo["freshRegions"] = [SimpleRegion(pl.coordPixels, (-1, -1, 2, 2)) for pl in toAdd]

    def SetFitPixelsRange(self, minFitPixels, maxFitPixels):
        self.logger.info("SetFitPixelsRange: %d %d" % (minFitPixels, maxFitPixels))
        self._minFitPixels = minFitPixels
        self._maxFitPixels = maxFitPixels
        
    def GetBestLocFit(self):
        self.logger.info("GetBestLocFit")
        res = None
        for potentialLoc in self._potentialLocs:
            locFitTmp = potentialLoc.GetBestLocFit()
            if locFitTmp is None:
                continue
            
            if res is None or locFitTmp.goodnessValue > res.goodnessValue:
                res = locFitTmp
        return res
    
    def AddLocationAndUpdate(self, locFit):
        self.logger.info("AddLocationAndUpdate")
        
        # Add to results and remove from image
        self.resultLocs.append(locFit)
        self._frame.SubtractLocs(self._psf, [locFit])
        
        # Save iteration info
        if self.saveIterationInfo:
            iterationInfo = self._localizer._iterationInfo
            iterationInfo["frameNext"] = self._frame.CopyFrame()
            iterationInfo["psfsToDelete"] = iterationInfo["startFrame"].data - iterationInfo["frameNext"].data
            iterationInfo["locs"] = [locFit]
            iterationInfo["oldRegions"] = [SimpleRegion(pl.coordPixels, (-1, -1, 2, 2)) for pl in self._potentialLocs]
        
        # Discard all potential locs that may be affected
        affectedBox = self._frame.GetBoxAround(locFit.coordPixels, \
            (locFit.bboxHWx + self._maxFitPixels, locFit.bboxHWy + self._maxFitPixels))
        self._potentialLocs = [potentialLoc for potentialLoc in self._potentialLocs if not potentialLoc.IsAffected(affectedBox)]
        
        # Analyze and add points from affected region
        self.AddPotentialLocs(bbox = affectedBox)
        
    @property
    def minFitPixels(self):
        return self._minFitPixels

    @property
    def maxFitPixels(self):
        return self._maxFitPixels

#===============================================================================
# PotentialLoc
#===============================================================================

class PotentialLoc(object):
    
    def __init__(self, collection, coord):
        self.logger = collection._localizer.logger
        self._collection = collection
        self.coord = coord
        self.coordPixels = collection._frame.GetCoordsPixels(self.coord)
        
        # Calc
        xI, yI = self.coordPixels
        directlyAroundData = collection._frame.data[(xI - 1):(xI + 2), (yI - 1):(yI + 2)]
        
        # Area
        self.area = (directlyAroundData > self._collection._localizer._noiselevel).sum()
        
        # Centroid
        tmpX, tmpY = np.meshgrid(np.arange(directlyAroundData.shape[0]), \
                         np.arange(directlyAroundData.shape[1]), \
                         indexing = "ij")
        centroidXI = np.average(tmpX, weights = directlyAroundData) - 1.0 + xI
        centroidYI = np.average(tmpY, weights = directlyAroundData) - 1.0 + yI
        self.weightedCentroid = collection._frame.GetPixelCoords((centroidXI, centroidYI))
            
        self._fitCache = {}
        
    def GetBestLocFit(self):
        res = None
        for fitPixels in range(self._collection._minFitPixels, self._collection._maxFitPixels + 1):
            fitLoc = self._Fit(fitPixels)
            if fitLoc is None:
                continue
            
            if res is None or fitLoc.goodnessValue > res.goodnessValue:
                res = fitLoc
        return res
    
    def IsAffected(self, bbox):
        if self.coordPixels[0] < bbox[0] or self.coordPixels[0] >= bbox[2]:
            return False
        if self.coordPixels[1] < bbox[1] or self.coordPixels[1] >= bbox[3]:
            return False
        return True
    
    def _Fit(self, fitPixels):
        
        def CheckForDuplicates(locFit):
            if locFit is None:
                return None
            
            for loc in self._collection.resultLocs:
                distXY = np.sqrt((loc.x - locFit.x) ** 2.0 + (loc.y - locFit.y) ** 2.0)
                if distXY <= localizer._duplicateMinDist:
                    self.logger.debug("Dismissed, duplicate dist %.2f nm" % (1e9 * distXY))
                    return None
            return locFit
        
        # Init
        localizer = self._collection._localizer
        frame = self._collection._frame
        psf = self._collection._psf
        axialCalibrator = self._collection._axialCalibrator
        fitMode = self._collection._fitMode
        
        # Check cache
        if fitPixels in self._fitCache:
            locFit = CheckForDuplicates(self._fitCache[fitPixels])
            self._fitCache[fitPixels] = locFit
            return locFit        
        self._fitCache[fitPixels] = None

        # Pick initial fit coord
        initialCoord = None
        if localizer._initialFitLoc == "localMaxima":
            initialCoord = self.coord
        elif localizer._initialFitLoc == "weightedCentroid":
            initialCoord = self.weightedCentroid
        else:
            raise NotImplementedError()
             
        # Do fitting
        locFit = psf.Fit(frame, initialCoord, fitPixels, \
                         fitMode = fitMode, \
                         initialMaxDist = localizer._initialMaxDist)
        
        if locFit is None:
            return None
                
        if min(locFit.x, locFit.y) < localizer._borderRegionWidth or \
            locFit.x + localizer._borderRegionWidth > frame.sizeX or \
            locFit.y + localizer._borderRegionWidth > frame.sizeY:
            self.logger.debug("Dismissed, in border requion")
            return None
        
        if locFit.z is not None:
            if locFit.z < axialCalibrator._calibFrom or locFit.z > axialCalibrator._calibTo:
                self.logger.debug("Dismissed, outside axial calibration range %.2f nm" % \
                                  (1e9 * locFit.z))
                return None

        if locFit.goodnessValue < localizer._minGoodness:
            self.logger.debug("Dismissed, to bad goodness %.2f" % (locFit.goodnessValue))
            return None
        
        if locFit.offset < localizer._minOffset:
            self.logger.debug("Dismissed, to bad offset %.2f" % (locFit.offset))
            return None
               
        locFit = CheckForDuplicates(locFit)
        
        # Save additional info
        if locFit is not None:
            locFit.area = self.area
            # TODO:
            locFit.distFitWeightedCentroid = np.sqrt((self.weightedCentroid[0] - locFit.bestFitParams[1]) ** 2.0 + \
                                                     (self.weightedCentroid[1] - locFit.bestFitParams[2]) ** 2.0 )
        
        # Cache
        self._fitCache[fitPixels] = locFit
        return locFit
        

if __name__ == "__main__":
    pass
