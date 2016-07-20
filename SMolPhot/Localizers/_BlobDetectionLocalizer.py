"""Implements blob detection localizer.

"""
import logging
import numpy as np
from SMolPhot.Components.BaseClasses import Param
from _Common import LocalizerBase
from skimage.feature import blob_log

__all__ = ["BlobDetectionLocalizer"]


class BlobDetectionLocalizer(LocalizerBase):
    """Blob detection localizer.

    Params:
        detThreshold (float): threshold = mean + detThreshold * std
        minSigma (float): minimum std for Gaussian kernel in px
        maxSigma (float): maximum std for Gaussian kernel in px
        fitPixels (int): half-width of the rectangular fit region

    """

    def __init__(self, params=[], **kwargs):
        self.logger = logging.getLogger("SMolPhot.Localizers.BlobDetectionLocalizer")
        paramsThis = [Param("maxIterations", 5,
                            friendlyName="Maximum # of iterations",
                            guiLimits=(1, 100)),

                      Param("detThreshold", 200.0,
                            friendlyName="Maxima threshold",
                            guiSuffix=" ph",
                            guiStep=1.0,
                            guiLimits=(1.0, 100000)),

                      Param("minSigma", 1.0,
                            friendlyName="Min sigma",
                            guiStep=0.1,
                            guiSuffix=" px"),

                      Param("maxSigma", 5.0,
                            friendlyName="Max sigma",
                            guiStep=0.1,
                            guiSuffix=" px"),

                      Param("numSigma", 10,
                            friendlyName="# of sigmas",
                            guiStep=1),

                      Param("fitPixels", 1,
                            friendlyName="Region half - width",
                            guiSuffix=" px",
                            guiLimits=(1, 10)),

                      Param("_initialMaxDist", 200e-9,
                            friendlyName="Initial max dist",
                            guiSiPrefix=True, guiSuffix="m",
                            guiStep=10e-9),

                      Param("_borderRegionWidth",
                            300e-9,
                            friendlyName="Border region",
                            guiSiPrefix=True,
                            guiSuffix="m",
                            guiStep=10e-9,
                            guiLimits=(0.0, 1.0)),

                      Param("_minPhotons",
                            1000.0,
                            friendlyName="Min photons",
                            guiSuffix=" ph",
                            guiDecimals=4,
                            guiStep=10.0,
                            guiLimits=(1.0, 1e6)),

                      Param("_duplicateMinDist",
                            200e-9,
                            friendlyName="Duplicate min dist",
                            guiSiPrefix=True,
                            guiSuffix="m",
                            guiStep=10.0,
                            guiLimits=(0.0, 1000e-9)),

                      Param("_removeOffset", 0,
                            friendlyName="Remove offset")]  # for testing only (0 or 1)
        LocalizerBase.__init__(self, params + paramsThis, **kwargs)

    def OptimizeForSeries(self, fseries):
        pass

    def FindMolecules(self, frame, psf, axialCalibrator, fitMode="z"):
        self.logger.info("Find molecules frame #%d, fitmode %s" % (frame.nr, fitMode))
        frameData = frame.data.copy()
        locs = []

        for i in range(self.maxIterations):
            locsAdded = []
            blobs = blob_log(frameData, min_sigma=self.minSigma, max_sigma=self.maxSigma, num_sigma=self.numSigma,
                             threshold=self.detThreshold)

            if len(blobs) == 0: # No molecules were found
                break

            maxIndices = blobs[:,:2].T.astype('uint8')
            maxCoords = frame.GetPixelCoords(maxIndices)

            self.logger.debug("Max Indices %s" % (str(maxIndices)))
            #self.logger.debug("Max Indices values %s" % (str(frame.data[maxIndices]))) TODO non integer indices

            pixelsX, pixelsY = frame.pixelsX, frame.pixelsY

            for x, y, xI, yI in zip(maxCoords[0], maxCoords[1], maxIndices[0], maxIndices[1]):
                self.logger.debug("Potential local maximum at pixel (%d, %d)" % (xI, yI))
                if min(xI, yI) < 1 or \
                    xI + 1 >= pixelsX or \
                    yI + 1 >= pixelsY:
                    self.logger.debug("Dismissed, pixel on edge")
                    continue

                locFit = psf.Fit(frame, (x, y), self.fitPixels, \
                             fitMode = fitMode, \
                             initialMaxDist = self._initialMaxDist)

                if locFit is None:
                    continue

                if min(locFit.x, locFit.y) < self._borderRegionWidth or \
                    locFit.x + self._borderRegionWidth > frame.sizeX or \
                    locFit.y + self._borderRegionWidth > frame.sizeY:
                    self.logger.debug("Dismissed, in border requion")
                    continue

                if np.sqrt((locFit.x - x) ** 2.0 + (locFit.y - y) ** 2.0) > self._initialMaxDist:
                    self.logger.debug("Dismissed, far from origin")
                    continue

                if locFit.photons < self._minPhotons:
                    self.logger.debug("Dismissed, photon flux low %.2f" % (locFit.photons))
                    continue

                if locFit.z is not None:
                    if locFit.z < axialCalibrator._calibFrom or locFit.z > axialCalibrator._calibTo:
                        self.logger.debug("Dismissed, outside axial calibration range %.2f nm" % (1e9 * locFit.z))
                        continue

                locsAdded.append(locFit)

            if len(locsAdded) == 0: # No molecules were found
                break
            locs += locsAdded


            # Check for duplicates
            duplicate = [False] * len(locs)
            for i in range(len(locs)):
                for j in range(i + 1, len(locs)):
                    if duplicate[j] or duplicate[i]:
                        continue
                    distXY = np.sqrt((locs[i].x - locs[j].x) ** 2.0 + (locs[i].y - locs[j].y) ** 2.0)
                    if distXY < self._duplicateMinDist:
                        # More then one molecule at same locationd
                        duplicate[j] = True
                        self.logger.debug("Found duplicate locations, removing one.")

            for i in range(len(locs) - 1, -1, -1):
                if duplicate[i]:
                    del locs[i]
                    self.logger.debug("Deleted loc %d" % (i))

            # Residual
            for loc in locsAdded:
                xs, ys = frame.coords
                fitParams = list(loc.bestFitParams)
                if self._removeOffset == 1:
                    fitParams[-1] = 0  # Assumes that offset is the last parameter and sets it to 0
                psfData = psf.Calc((xs, ys), *fitParams)
                frameData -= psfData

        return locs

    @property
    def name(self):
        return "Blob Detection Localizer"


if __name__ == '__main__':
    pass