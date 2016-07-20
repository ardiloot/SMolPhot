"""Implements Frame and FrameSeries classes for loading frames.

"""

import logging
import yaml
import numpy as np
from os import path
from PIL import Image
from BaseClasses import ParamsBaseClass, Param
from collections import defaultdict

__all__ = ["Frame", "FrameSeries"]

#===============================================================================
# FrameSeries
#===============================================================================

class FrameSeries(ParamsBaseClass):
    """Class for a series of frames. Usually created by static method 
    FromMetafile(), which loads all required parameters from special yaml file.
        
    Args:
        pixelSize (float): single pixel width (square pixels are assumed)
        frameFiles (list of str): the list of the files of frames
        stackFile (str): the filename of the TIFF stack  
        params (list of Params): unused, needed if inherited
        kwargs: used to set the values of the parameters
        
    """
        
    def __init__(self, pixelSize, series = "", frameFiles = None, stackFile = None, \
                 params = [], maxFramesToLoad = None, \
                 activationsFile = None, preprocessors = [], \
                 excludedBorderWidth = 300e-9, **kwargs):
        
        self.logger = logging.getLogger("SMolPhot.FrameSeries")
        self.logger.debug("%s: __init__" % (series))
        
        self.pixelSize = pixelSize
        self._series = series
        self._frameFiles = frameFiles
        self._stackFile = stackFile
        self._maxFramesToLoad = int(1e9) if maxFramesToLoad is None else maxFramesToLoad
        self._activationsFile = activationsFile
        self._preprocessors = preprocessors
        self._excludedBorderWidth = excludedBorderWidth
        
        self._metadatafile = ""
        self._origninalFrames = []
        self._preprocessedFrames = []
        self._zs = None
        self._precalcPixelShape = (0, 0)
        self._precalcCoords = None
        

        if self._frameFiles is not None and self._stackFile is not None:
            raise ValueError("Only one, frameFiles or stackFile must be given.")
        
        paramsThis = [Param("name", "test dataset"),
                      Param("cameraQE", None, required = False),
                      Param("wl", 660e-9),
                      Param("NA", None, required = False),
                      Param("readoutGaussianNoise", None, required = False),
                      Param("emGainGamma", None, required = False),
                      Param("spuriousNoisePoisson", None, required = False),
                      Param("totalGain", None),
                      Param("adcElectronConversion", None, required = False),
                      Param("adcBaseline", None),
                      Param("adcSaturation", None, required = False),
                      Param("adcQuantization", None, required = False),
                      Param("zRangeMin", None, required = False),
                      Param("zRangeMax", None, required = False),
                      Param("zStep", None, required = False)]
        
        ParamsBaseClass.__init__(self, params + paramsThis, **kwargs)
        
        # Check if, we have all we need
        self.AssertRequiredParamsSet()
        
        self._LoadFrames()
        self._LoadGroundTruth()
        
        self.logger.debug("%s: __init__ done" % (self._series))
        
         
    def _LoadFrames(self):
        
        def ConvertADU2Photons(data):
            res = data - self.adcBaseline
            res /= self.totalGain
            return res
        
        self.logger.info("LoadFrames from series '%s'" % (self._series))
        self._origninalFrames = []
        self._preprocessedFrames = []
        
        if self._frameFiles is not None:
            # From list of files
            for frameNr in len(self._frameFiles):
                filename = self._frameFiles[frameNr]
                data = ConvertADU2Photons(np.array(Image.open(filename), dtype = np.float).T)
                self._origninalFrames.append(Frame(self, frameNr, data))
                
        elif self._stackFile is not None:
            # Froms TIFF stack
            self._stackPIL = Image.open(self._stackFile)
            frameNr = 0
            try:
                while frameNr < self._maxFramesToLoad:
                    self._stackPIL.seek(frameNr)
                    data = ConvertADU2Photons(np.array(self._stackPIL, dtype = np.float).T)
                    self._origninalFrames.append(Frame(self, frameNr, data))
                    frameNr += 1
            except EOFError:
                pass
        else:
            raise ValueError("frameFiles or stackFile must be given.")
            
            
        # Progressed frames
        self._preprocessedFrames = [None] * len(self)
            
        # Update zs-values
        if self.zRangeMin is not None:
            self._zs = np.arange(self.zRangeMin, self.zRangeMax + 1e-14, self.zStep)
        
        self.logger.info("LoadFrames series '%s' done (frame count %d)" % (self._series, len(self)))
        
    def _LoadGroundTruth(self):    
        self.logger.info("LoadGroundTruth series '%s'" % (self._series))
        
        # Activations
        if self._activationsFile is not None:
            #profile = cProfile.Profile()
            #profile.enable()
            
            nrs, frameNrs, xs, ys, zs, indensities = \
                np.loadtxt(self._activationsFile, \
                           delimiter = ",", \
                           unpack = True, \
                           skiprows = 1)
            
            # Discard ground-truth near the border
            groundTruthToKeep = []
            frame = self.GetOriginalFrame(0)
            sizeX, sizeY = frame.sizeX, frame.sizeY
            for i in range(len(frameNrs)):
                x, y = 1e-9 * xs[i], 1e-9 * ys[i]
                if min(x, y) < self._excludedBorderWidth or \
                    x + self._excludedBorderWidth > sizeX or \
                    y + self._excludedBorderWidth > sizeY:
                    continue
                groundTruthToKeep.append(i)
            groundTruthToKeep = np.array(groundTruthToKeep)
            
            # Save ground-truth variables
            self._actNr = np.ascontiguousarray(nrs[groundTruthToKeep].astype(int))
            self._actFrameNr = np.ascontiguousarray(frameNrs[groundTruthToKeep].astype(int) - 1)
            self._actCoords = np.ascontiguousarray(1e-9 * xs[groundTruthToKeep]), \
                              np.ascontiguousarray(1e-9 * ys[groundTruthToKeep]), \
                              np.ascontiguousarray(1e-9 * zs[groundTruthToKeep])
            self._actCoordsArray = np.ascontiguousarray(np.array(self._actCoords).T)
            self._actIs = np.ascontiguousarray(indensities[groundTruthToKeep])
            
            discardedPoints = len(nrs) - len(self._actNr)
            if discardedPoints > 0:
                self.logger.warn("Discarded %d ground-truth points, because of excluded boundary width %.1f nm" % \
                                 (discardedPoints, 1e9 * self._excludedBorderWidth))
                
            # Build frame map
            self._actFrameMap = defaultdict(list)
            for i in np.arange(len(self._actFrameNr)):
                self._actFrameMap[self._actFrameNr[i]].append(i)
            
            # Convert to numpy arrays
            for k, v in self._actFrameMap.iteritems():
                self._actFrameMap[k] = np.ascontiguousarray(np.array(v, dtype = int))
    
            # Number molecules
            dictActMolIdsTmp = dict()
            curMolIdTmp = 0
            self._actMolId = np.zeros_like(self._actFrameNr, dtype = int)
            
            for i in np.arange(len(self._actFrameNr)):
                coord = (xs[i], ys[i], zs[i])
                if coord in dictActMolIdsTmp:
                    self._actMolId[i] = dictActMolIdsTmp[coord]
                else:
                    self._actMolId[i] = curMolIdTmp
                    dictActMolIdsTmp[coord] = curMolIdTmp
                    curMolIdTmp += 1

            #profile.disable()
            #pstats.Stats(profile).sort_stats("cumulative").print_stats(15)

        self.logger.info("LoadGroundTruth done series '%s'" % (self._series))
        
    @staticmethod
    def FromMetafile(filename, preprocessors, series = "sequence", maxFrameOverRide = None):
        # series: sequence or axial calibration
        with open(filename, "r") as stream:
            metadata = yaml.safe_load(stream)

        if not series in metadata["series"]:
            raise ValueError("Could not find series %s" % (series))

        dataDir = path.dirname(path.abspath(filename))
        pixelSize = metadata["pixel size"]
        paramValues = metadata["frame params"]
        seriesMetadata = metadata["series"][series]
        
        if "max frames" in seriesMetadata and maxFrameOverRide is None:
            paramValues["maxFramesToLoad"] = seriesMetadata["max frames"]
        elif maxFrameOverRide is not None:
            paramValues["maxFramesToLoad"] = maxFrameOverRide

        if "ground-truth" in seriesMetadata:
            groundTruthMetadata = seriesMetadata["ground-truth"]
            activationsFile = path.join(dataDir, groundTruthMetadata["activations"])         
            paramValues["activationsFile"] = activationsFile
            
            if "excludedBorderWidth" in groundTruthMetadata:
                paramValues["excludedBorderWidth"] = groundTruthMetadata["excludedBorderWidth"]

        if "params" in seriesMetadata:
            paramValues.update(seriesMetadata["params"])

        # Init FrameSeries
        if "qfiles" in seriesMetadata:
            frameFiles = [path.join(dataDir, f) for f in seriesMetadata["qfiles"]]
            res = FrameSeries(pixelSize, series = series, frameFiles = frameFiles, \
                              preprocessors = preprocessors, **paramValues)
        elif "stack file" in seriesMetadata:
            stackFile = path.join(dataDir, seriesMetadata["stack file"])
            res = FrameSeries(pixelSize, series = series, stackFile = stackFile, \
                              preprocessors = preprocessors, **paramValues)
        else:
            raise ValueError("Frame files or stackFile must be given.")
        
        res._metadatafile = filename
        return res
    
    def GetPixelCoords(self, indices):
        indicesX, indicesY = indices
        # Center coordinates
        coordX = self.pixelSize * (0.5 + indicesX)
        coordY = self.pixelSize * (0.5 + indicesY)
        return coordX, coordY
    
    def GetCoordsPixels(self, coords):
        if isinstance(coords[0], np.ndarray) and len(coords[0]) > 1:
            indicesX = np.maximum(0, np.around(coords[0] / self.pixelSize - 0.5).astype(int))
            indicesY = np.maximum(0, np.around(coords[1] / self.pixelSize - 0.5).astype(int))
        else:
            indicesX = max(0, int(round(float(coords[0]) / self.pixelSize - 0.5)))
            indicesY = max(0, int(round(float(coords[1]) / self.pixelSize - 0.5)))
        return indicesX, indicesY
    
    def GetCoordsMatrix(self, pixelShape):
        if self._precalcPixelShape != pixelShape:
            self.logger.debug("GetCoordsMatrix(%s) updated" % (str(pixelShape)))
            pixelsX, pixelsY = np.meshgrid(np.arange(pixelShape[0]), \
                                     np.arange(pixelShape[1]), indexing = "ij")
            self._precalcCoords = self.GetPixelCoords((pixelsX, pixelsY))
            self._precalcPixelShape = pixelShape
        return self._precalcCoords
    
    def SaveMetadata(self, filename):
        metadata = {}
        metadataParams = metadata["frame params"] = {} 

        for param, value in self.GetParams(asList = True):
            metadataParams[param.name] = value
            
        metadata["frame pixel size"] = self.pixelSize            
        metadata["frame qfiles"] = [path.basename(f) for f in self._frameFiles]    
        
        with open(filename, "w") as stream:
            yaml.safe_dump(metadata, stream)
    
    def GetOriginalFrame(self, nr):
        self.logger.debug("%s: GetOriginalFrame(%d)" % (self._series, nr))
        return self._origninalFrames[nr]
    
    def OptimizePreprocessors(self):
        self.logger.debug("%s: OptimizePreprocessors" % (self._series))
        for preprocessor in self._preprocessors:
            preprocessor.OptimizeForSeries(self)
        self.logger.debug("%s: OptimizePreprocessors done" % (self._series))
        
    def GetPreprocessedFrame(self, nr, roi = None):
        # TODO: performance issue, avoid recalculation if possible
        self.logger.debug("%s: GetPreprocessedFrame(%d)" % (self._series, nr))
        newFrame = self._origninalFrames[nr].CopyFrame()
        
        for preprocessor in self._preprocessors:
            preprocessor.Apply(newFrame)
        
        newFrame.Crop(roi)
        return newFrame
    
    def HasGroundTruth(self):
        return self._activationsFile is not None
    
    def GetNumberOfActMolecules(self, frameFrom = 0, frameTo = None):
        if self is None or not self.HasGroundTruth():
            return None
        
        if frameTo is None:
            frameTo = len(self)
        
        res = 0
        for frameNr in range(frameFrom, frameTo):
            res += len(self._actFrameMap[frameNr])
        return res

            
    def __len__(self):
            return len(self._origninalFrames)
    
    @property
    def series(self):
        return self._series
    
    @property
    def metadatafile(self):
        return self._metadatafile

    @property
    def activationsFile(self):
        return self._activationsFile
    
    @property
    def zs(self):
        return self._zs        


#===============================================================================
# Frame
#===============================================================================

class Frame(ParamsBaseClass):
    """Class for a single frame. Usually creaded by FrameSeries or copyed by
    CopyFrame() method. Has some precalc features.
        
    Args:
        fseries (FrameSeries): reference to frame series 
        nr (int): number of frame in fseries (not frame number form filename!)
        data (numpy array): frame data
        
    """
    
    def __init__(self, fseries, nr, data):
        self.logger = logging.getLogger("SMolPhot.Frame")
        self.logger.debug("%s-%d: __init__" % (fseries.series, nr))
        
        self._fseries = fseries
        self._nr = nr
        self.data = data
        self.logger.debug("%s-%d: __init__ done" % (fseries.series, nr))

    # Public methods
       
    def CopyFrame(self, roi = None):
        if roi is None:
            data = self.data
        else:
            (x0, x1), (y0, y1) = roi
            data = self.data[x0:x1, y0:y1]
            
        res = Frame(self._fseries, self.nr, data.copy())
        return res
    
    def Crop(self, roi):
        if roi == None:
            return
        
        (x0, x1), (y0, y1) = roi
        self.data = self.data[x0:x1, y0:y1]
       
    def GetPixelCoords(self, *args, **kwargs):
        return self._fseries.GetPixelCoords(*args, **kwargs)
    
    def GetCoordsPixels(self, *args, **kwargs):
        return self._fseries.GetCoordsPixels(*args, **kwargs)
   
    def HasGroundTruth(self):
        return self._fseries.HasGroundTruth()
   
    def GetActMolecules(self):
        indices = self._fseries._actFrameMap[self.nr]
        cX, cY, cZ = self._fseries._actCoords
        return cX[indices], cY[indices], cZ[indices], self._fseries._actIs[indices]

    def IsPixelInside(self, (xI, yI)):
        if xI < 0 or yI < 0 or xI >= self.data.shape[0] or yI >= self.data.shape[1]:
            return False
        else:
            return True
        
    def CropCoord(self, (xI, yI)):
        res = (max(0, min(self.data.shape[0] - 1, xI)), max(0, min(self.data.shape[1] - 1, yI)))
        return res
        
    def FindMaximas(self, threshold, bbox = None, localMaxima = True):
        pixelsX, pixelsY = self.pixelsX, self.pixelsY
        
        if bbox is None:
            sl = (slice(0, pixelsX), slice(0, pixelsY)) 
        else:
            sl = (slice(bbox[0], bbox[2]), slice(bbox[1], bbox[3]))
                
        maxIndicesInitialX, maxIndicesInitialY  = np.where(self.data[sl] > threshold)
        if bbox is not None:
            maxIndicesInitialX += bbox[0]
            maxIndicesInitialY += bbox[1]
        
        res = []
        for xI, yI in zip(maxIndicesInitialX, maxIndicesInitialY):
            if min(xI, yI) < 1 or xI + 1 >= pixelsX or yI + 1 >= pixelsY:
                #self.logger.debug("Dismissed, pixel on edge")
                continue
            
            if localMaxima:
                # Local maxima is required
                pixelValue = self._data[xI, yI]
                directlyAroundData = self._data[(xI - 1):(xI + 2), (yI - 1):(yI + 2)]
                maxAround = np.max(directlyAroundData)
                
                if maxAround > pixelValue:
                    # Not local maxima
                    continue
                
            coord = self.GetPixelCoords((xI, yI))
            res.append(coord)
                
        return res
    
    def GetBoxAround(self, (xI, yI), (hwX, hwY)):
        #print type(xI), type(hwX)
        x0 = max(0, int(round(xI - hwX)))
        x1 = min(int(round(xI + hwX + 1)), self.data.shape[0])
        y0 = max(0, int(round(yI - hwY)))
        y1 = min(int(round(yI + hwY + 1)), self.data.shape[1])
        return x0, y0, x1, y1

    def GetAround(self, (xI, yI), halfPixelsAround, flatten = False):
        x0, y0, x1, y1 = self.GetBoxAround((xI, yI), (halfPixelsAround, halfPixelsAround))

        coordsAround = (self.coords[0][x0:x1, y0:y1], self.coords[1][x0:x1, y0:y1])
        dataAround = self._data[x0:x1, y0:y1]
        molCoordAround = (xI - x0, yI - y0)
        
        if flatten:
            return molCoordAround, (coordsAround[0].ravel(), coordsAround[1].ravel()), \
                dataAround.ravel()
        else:
            return molCoordAround, coordsAround, dataAround

    def SubtractLocs(self, psf, locs):
        for loc in locs:
            # TODO: subraction only near the molecule
            self.data -= psf.CalcWithoutOffset(self.coords, *loc.bestFitParams) 
        return self
    
    # Getters and setters
       
    @property
    def data(self):
        return self._data
        
    @data.setter
    def data(self, value):
        self.logger.debug("set data #%d" % (self._nr))
        self._mean = None
        self._std = None
        self._data = value
    
    @property
    def filename(self):
        self.logger.debug("get filename #%d" % (self._nr))
        return self._filename
    
    @property    
    def nr(self):
        return self._nr
    
    @property
    def mean(self):
        self.logger.debug("get mean #%d" % (self._nr))
        if self._mean is None:
            self._mean = np.mean(self._data)
        return self._mean
    
    @property
    def std(self):
        self.logger.debug("get std #%d" % (self._nr))
        if self._std is None:
            self._std = np.std(self._data)
        return self._std

    @property
    def pixelsX(self):
        return self.data.shape[0]    
    
    @property
    def pixelsY(self):
        return self.data.shape[1]

    @property
    def sizeX(self):
        return self.pixelsX * self._fseries.pixelSize    
    
    @property
    def sizeY(self):
        return self.pixelsY * self._fseries.pixelSize
    
    @property
    def z(self):
        return self._fseries._zs[self.nr]
    
    @property
    def coords(self):
        return self._fseries.GetCoordsMatrix(self.data.shape)

    @property
    def pixelSize(self):
        return self._fseries.pixelSize

if __name__ == "__main__":
    pass