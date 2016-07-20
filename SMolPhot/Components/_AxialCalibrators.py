import logging
import numpy as np
from BaseClasses import ParamsBaseClass, Param
# import cProfile, pstats

__all__ = ["ManualAxialCalibrator",
           "GroundTruthAxialCalibrator"]

#===============================================================================
# AxialCalibratorBase
#===============================================================================

class AxialCalibratorBase(ParamsBaseClass):
    
    def __init__(self, params = [], **kwargs):
        self._fittedLocs = None
        
        paramsThis = [Param("_calibFrom", -1000e-9,
                            friendlyName = "Calib from z",
                            guiSiPrefix = True,
                            guiSuffix = "m",
                            guiStep = 10e-9,
                            guiLimits = (-10e-6, 10e-6)),
                      
                      Param("_calibTo", 1000e-9,
                            friendlyName = "Calib to z",
                            guiSiPrefix = True,
                            guiSuffix = "m",
                            guiStep = 10e-9,
                            guiLimits = (-10e-6, 10e-6)),
                      ]
        
        ParamsBaseClass.__init__(self, params + paramsThis, **kwargs)
    
    @property
    def name(self):
        raise NotImplementedError()

#===============================================================================
# ManualAxialCalibrator
#===============================================================================

class ManualAxialCalibrator(AxialCalibratorBase):
    """Simple axial calibrator fitts PSF to every centerpoint defined by ROI.
    If the fitted location is to far away, it discards it. Currently only
    information form the first ROI is used.
    
    """
    
    def __init__(self, params = [], **kwargs):
        self.logger = logging.getLogger("SMolPhot.AxialCalibrators.ManualAxialCalibrator")
        
        paramsThis = [Param("_minFitPixels", 3,
                            friendlyName = "Min fit half-width",
                            guiSuffix = " px",
                            guiLimits = (1, 10)),
                      
                      Param("_maxFitPixels", 5,
                            friendlyName = "Max fit half-width",
                            guiSuffix = " px",
                            guiLimits = (1, 10)),
                      
                      Param("_initialMaxDist", 300e-9,
                            friendlyName = "Initial max dist",
                            guiSiPrefix = True,
                            guiSuffix = "m",
                            guiStep = 10e-9),
                      ]
        
        AxialCalibratorBase.__init__(self, params + paramsThis, **kwargs)
        
    def CalibratePsf(self, psf, axialFseries, rois):
        """Calibrates PSF
        
        Args:
            psf: psf to calibrate (is also used for fitting)
            axialFseries (FrameSeries): axial calibration frame series
            rois (list of MyEllipseROI): list of ROIs
        """
        
        self.logger.debug("CalibratePsf: ROIS: %s" % (str(rois)))
        
        # Exctract center points and angles from elliptical rois
        initPosXs = np.array([state["pos"][0] + 0.5 * state["size"][0] for state in rois])
        initPosYs = np.array([state["pos"][1] + 0.5 * state["size"][1] for state in rois])
        # initAngles = np.array([state["angle"] for state in rois])
        
        startPositions = [zip(initPosXs, initPosYs) for _ in range(len(axialFseries))]
        self._DoCalibrationFromLocs(psf, axialFseries, startPositions)
        
    def _DoCalibrationFromLocs(self, psf, axialFseries, startPositions):
        # Start profiling
        # profile = cProfile.Profile()
        # profile.enable()
      
        # For every frame fit every molecule
        calibrationPoints = []
        for frameNr in range(len(axialFseries)):
            self.logger.debug("Processing frame %d" % (frameNr))
            
            # Load frame
            frame = axialFseries.GetPreprocessedFrame(frameNr)
            
            if frame.z < self._calibFrom or frame.z > self._calibTo:
                self.logger.debug("Outside calibration reqion")
                continue
            
            # Fit every molecule
            fitResThisFrame = []
            for i, (startX, startY) in enumerate(startPositions[frameNr]):
                self.logger.debug("Processing molecule %d" % (i))
                
                # Vary fitregion size from minFitPixels to maxFitPixels
                for fitPixels in range(self._minFitPixels, self._maxFitPixels + 1):
                    locFit = psf.Fit(frame, (startX, startY), fitPixels, \
                                     fitMode = "sigma", \
                                     initialMaxDist = self._initialMaxDist)
                    
                    if locFit is None:
                        self.logger.warn("No location was found")
                        continue
                    
                    locFit.startPositionNr = i
                    locFit.startPositionOffset = (startX - locFit.x, startY - locFit.y)
                    
                    fitResThisFrame.append(locFit)
                    
            # Add fitter results from this frame to calibration points
            calibrationPoints.append((frame, fitResThisFrame))
       
        psf.SetAxialCalibPoints(calibrationPoints, useStds = True)
                
        # profile.disable()
        # pstats.Stats(profile).sort_stats("cumulative").print_stats(15) 
        
    @property
    def name(self):
        return "Manual calibrator"
        

#===============================================================================
# GroundTruthAxialCalibrator
#===============================================================================

class GroundTruthAxialCalibrator(ManualAxialCalibrator):
    """Interpolating axial calibrator fits all molecules to produce extra
    calibration points
    
    """
    
    def __init__(self, params = [], **kwargs):
        paramsThis = []
        ManualAxialCalibrator.__init__(self, params + paramsThis, **kwargs)
        self.logger = logging.getLogger("SMolPhot.AxialCalibrators.GroundTruthAxialCalibrator")
        
    def CalibratePsf(self, psf, axialFseries, rois):
        self.logger.info("CalibratePsf")
        
        startPositions = {}
        for frameNr in range(len(axialFseries)):
            gtPosX, gtPosY, _, _ = axialFseries.GetOriginalFrame(frameNr).GetActMolecules()
            startPositions[frameNr] = zip(gtPosX, gtPosY)
        
        self._DoCalibrationFromLocs(psf, axialFseries, startPositions)
         
    @property
    def name(self):
        return "Ground-truth calibrator"


if __name__ == '__main__':
    pass
