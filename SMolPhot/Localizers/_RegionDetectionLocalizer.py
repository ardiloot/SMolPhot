"""Implements region detection localizer.

"""
import logging
import numpy as np
from SMolPhot.Components.BaseClasses import Param
from _Common import LocalizerBase, CalcMinimumDistance, AddLocsWoDuplicates
from skimage.measure import label, regionprops

from collections import OrderedDict

__all__ = ["RegionDetectionLocalizer"]


def divideRegion(region, maxArea):
    if region.area <= maxArea:
        offset = np.array([0, 0])  # coordinate offset in pixels
        return [region], [offset]
    else:
        threshold = region.mean_intensity + 0.01 * (region.max_intensity - region.mean_intensity)
        mask = region.intensity_image > threshold
        labelImage = label(mask)
        subregions = regionprops(labelImage, intensity_image = region.intensity_image)
        regionList = []
        offsetList = []
        for subregion in subregions:
            newSubregions, newOffsets = divideRegion(subregion, maxArea)
            for newOffset in newOffsets:
                newOffset += np.array(region.bbox[:2])
            regionList += newSubregions
            offsetList += newOffsets
        return regionList, offsetList


def findHighestSubregion(region, maxArea):
    if region.area <= maxArea:
        offset = np.array([0, 0])  # coordinate offset in pixels
        return region, offset
    else:
        threshold = (region.max_intensity + region.mean_intensity) / 2.0
        mask = region.intensity_image > threshold
        labelImage = label(mask)
        subregions = regionprops(labelImage, intensity_image = region.intensity_image)
        maxI = region.min_intensity
        for subregion in subregions:
            if subregion.max_intensity > maxI:
                maxI = subregion.max_intensity
                highestSubregion = subregion
        newRegion, offset = findHighestSubregion(highestSubregion, maxArea)
        return newRegion, offset + np.array(region.bbox[:2])

class RegionDetectionLocalizer(LocalizerBase):

    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Localizers.RegionDetectionLocalizer")
        self._modeComboValues = OrderedDict([("Single peak", "Single peak"),
                                             ("Multipeak", "Multipeak")])

        paramsThis = [Param("_mode", "Single peak",
                            friendlyName = "Peak detection mode",
                            guiType = "list",
                            guiValues = self._modeComboValues.keys()),

                      Param("maxIterations", 5,
                            friendlyName = "Maximum # of iterations",
                            guiLimits = (1, 100)),

                      Param("detThreshold", 50.0,
                            friendlyName = "Maxima threshold",
                            guiSuffix = " ph",
                            guiStep = 1.0,
                            guiLimits = (1.0, 100000)),

                      Param("maskThreshold", 30.0,
                            friendlyName = "Labeling threshold",
                            guiSuffix = " ph",
                            guiStep = 1.0,
                            guiLimits = (1.0, 100000)),

                      Param("minArea", 4,
                            friendlyName = "Minimum area",
                            guiSuffix = " px"),

                      Param("maxArea", 30,
                            friendlyName = "Maximum area",
                            guiSuffix = " px"),

                      Param("fitPixels", 4,
                            friendlyName = "Region half-width",
                            guiSuffix = " px",
                            guiLimits = (1, 10)),
                      
                      Param("_minFitPixels", \
                            2, \
                            friendlyName = "Min region half-width", \
                            guiSuffix = " px", \
                            guiLimits = (1, 10)),

                      Param("_initialMaxDist", 100e-9,
                            friendlyName = "Initial max dist",
                            guiSiPrefix = True, guiSuffix = "m",
                            guiStep = 10e-9),

                      Param("_borderRegionWidth",
                            300e-9,
                            friendlyName = "Border region",
                            guiSiPrefix = True,
                            guiSuffix = "m",
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
                            guiStep = 0.1,
                            guiLimits = (0.0, 1e12))]
        

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
            
            return locFit
        
        def FindPotentialLocations(frame):
            mask = frame.data > self.maskThreshold
            labelImg = label(mask)
            regions = regionprops(labelImg, intensity_image = frame.data)
            internalRegions = []


            maxIndices = []
            if self._mode == "Single peak":
                for region in regions:
                    offset = np.array([0, 0])
                    region.offset = offset
                    
                    if region.area < self.minArea:
                        self.logger.debug("Region dismissed, small area (%d)" % (region.area))
                        continue
                    
                    if region.max_intensity < self.detThreshold:
                        self.logger.debug("Region dismissed, low maximum intensity (%f)" % (region.max_intensity))
                        continue
                
                    if region.area > self.maxArea:
                        region, offset = findHighestSubregion(region, self.maxArea)
                
                    maxIndices.append(np.array(region.weighted_centroid) + offset)
                    region.offset = offset
                    internalRegions.append(region)
                    # print "area", region.area, "convex area", region.convex_area, "orientation", region.orientation, "eccentricity", region.eccentricity, "indices", np.array(region.weighted_centroid) - offset

            elif self._mode == "Multipeak":
                for region in regions:
                    offset = np.array([0, 0])
                    region.offset = offset
                    
                    if region.area < self.minArea:
                        self.logger.debug("Region dismissed, small area (%d)" % (region.area))
                        continue
                    
                    if region.max_intensity < self.detThreshold:
                        self.logger.debug("Region dismissed, low maximum intensity (%f)" % (region.max_intensity))
                        continue
      
                    newRegions, newOffsets = None, None
                    if region.area > self.maxArea:
                        newRegions, newOffsets = divideRegion(region, self.maxArea)
                        for nr, no in zip(newRegions, newOffsets):
                            nr.offset = no
                            internalRegions.append(nr)
                    else:
                        newRegions, newOffsets = [region], [offset]
                        maxIndices.append(np.array(region.weighted_centroid) + offset)
                       
                    for newRegion, newOffset in zip(newRegions, newOffsets):
                        maxIndices.append(np.array(newRegion.weighted_centroid) + newOffset)
                       
            else:
                raise NotImplementedError()
            
            if len(maxIndices) == 0:
                maxIndices = np.zeros((0,)), np.zeros((0, ))
            else:
                maxIndices = np.array(maxIndices).T
                
            maxCoords = frame.GetPixelCoords(maxIndices)
            
            return maxIndices, maxCoords, [r for r in regions if r.area > self.minArea], internalRegions

        def RecursiveSearch(iteration, frame, fitPixels, locsAlreadyFound = []):
            if fitPixels < self._minFitPixels:
                return locsAlreadyFound
            
            iterationInfo = {"iteration": iteration, "fitPixels": fitPixels, \
                             "startFrame": frame}
            self._iterationsInfoTmp.append(iterationInfo)
            
            self.logger.info("RecursiveSearch iteration %d, fitPixels %d" % (iteration, fitPixels))

            _, maxCoords, regions, internalRegions = FindPotentialLocations(frame)
            iterationInfo["regions"] = regions
            iterationInfo["internalRegions"] = internalRegions
            
            locsThis = []
            for x, y in zip(maxCoords[0], maxCoords[1]):                
                locFit = _FitMolecule(frame, (x, y), fitPixels)
                
                if locFit is None:
                    continue
                
                locFit.iteration = iteration
                locsThis.append(locFit)
                
            locsThis.sort(key = lambda loc: -loc.goodnessValue)
            locs, locsThis = AddLocsWoDuplicates(locsAlreadyFound, locsThis, \
                                                 self._minDistSimFit,
                                                 self._duplicateMinDist, \
                                                 self.logger)                
            
            # Remove found PSFs
            frameNext = frame.CopyFrame().SubtractLocs(psf, locsThis)
            
            # Save additional iterationInfo
            iterationInfo["psfsToDelete"] = frame.data - frameNext.data
            iterationInfo["frameNext"] = frameNext
            iterationInfo["locs"] = locsThis
            
            fitPixelsNext = fitPixels - 1 if len(locsThis) == 0 else fitPixels
            
            self.logger.info("Locs from recursiveSearch iteration %d: %d" % (iteration, len(locsThis)))
            locs = RecursiveSearch(iteration + 1, frameNext, fitPixelsNext, locsAlreadyFound = locs)
            
            return locs
        
        #----------------------------------------------------------------------- 
        
        self._iterationsInfoTmp = []
        self.logger.info("Find molecules frame #%d, fitmode %s" % (inputFrame.nr, fitMode))
        locs = RecursiveSearch(0, inputFrame, self.fitPixels)
        locs = CalcMinimumDistance(locs)
        return locs

    @property
    def name(self):
        return "Region Detection Localizer"


if __name__ == '__main__':
    pass
