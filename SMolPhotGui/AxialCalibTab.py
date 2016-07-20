import Utils
import numpy as np
import pyqtgraph as pg
# import cProfile, pstats

from PyQt4 import QtGui
from collections import OrderedDict
from pyqtgraph.dockarea import Dock

#===============================================================================
# OriginalImageDock
#===============================================================================

class OriginalImageDock(Utils.ImageViewDock):
            
    def DrawPlot(self):
        # profile = cProfile.Profile()
        # profile.enable()
        
        axialFseries = self.main._axialFseries
        nr = self.tab._currentIndex
        
        # Load frame data
        originalFrame = axialFseries.GetOriginalFrame(nr)

        # Titles and labels
        self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
        self.label.setText("Orignial frame: #%d (z = %.3fum)" % \
                           (nr + 1, 1e6 * originalFrame.z))
        
        # Actual plotting
        transform = QtGui.QTransform()
        transform.scale(axialFseries.pixelSize, axialFseries.pixelSize)    
        self.imv.setImage(originalFrame.data, autoLevels = False, \
                          autoRange = False, transform = transform)

        # profile.disable()
        # pstats.Stats(profile).sort_stats("cumulative").print_stats(15) 

    def name(self):
        return "Original"

#===============================================================================
# PreprocessedImageDock
#===============================================================================

class PreprocessedImageDock(Utils.ImageViewDock):

    def DrawPlot(self):
        # profile = cProfile.Profile()
        # profile.enable()
        
        axialFseries = self.main._axialFseries
        nr = self.tab._currentIndex
        
        # Load frame data
        preprocessedFrame = axialFseries.GetPreprocessedFrame(nr)        
        
        # Ground truth
        if preprocessedFrame.HasGroundTruth():
            gtXs, gtYs, _, _ = preprocessedFrame.GetActMolecules()
            self.imv.plt.plot(gtXs, gtYs, symbol = "+", pen = None, \
                              symbolPen = "b", symbolBrush = None)
        
        # Titles and labels
        self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
        self.label.setText("Preprocessed frame: #%d (z = %.3fum)" % \
                           (nr + 1, 1e6 * preprocessedFrame.z))
     
        # Actual plotting
        transform = QtGui.QTransform()
        transform.scale(axialFseries.pixelSize, axialFseries.pixelSize)    
        self.imv.setImage(preprocessedFrame.data, autoLevels = False, \
                          autoRange = False, transform = transform)
        
        # profile.disable()
        # pstats.Stats(profile).sort_stats("cumulative").print_stats(15) 
        
    def name(self):
        return "Preprocessed"

#===============================================================================
# CalibrationDock
#===============================================================================

class CalibrationDock(Dock, Utils.SaveStateBaseClass):
    
    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = True
        
        Dock.__init__(self, *args, **kwargs)   
        
        self._axialPlot = pg.PlotWidget()
        self._axialPlot.addLegend()
        
        self._axialControls = QtGui.QWidget(self)
        self._axialControlsLayout = QtGui.QGridLayout(self._axialControls)
        
        self._axialButtons = QtGui.QWidget(self._axialControls)
        self._axialButtonsLayout = QtGui.QHBoxLayout(self._axialButtons)
        
        self._calibButton = QtGui.QPushButton("Calibrate", self._axialButtons)
        self._calibButton.setIcon(QtGui.QIcon(":/16x16/16x16/Play_16x16.png"))
        self._calibButton.pressed.connect(self._CalibrateAxialPosition)
        
        self._calibDelButton = QtGui.QPushButton("Clear", self._axialButtons)
        self._calibDelButton.setIcon(QtGui.QIcon(":/16x16/16x16/Delete_16x16.png"))
        self._calibDelButton.pressed.connect(self._DeleteAxialCalibration)
        
        self._axialButtonsLayout.addWidget(self._calibButton)
        self._axialButtonsLayout.addWidget(self._calibDelButton)
        
        self._axialControlsLayout.addWidget(self._axialPlot, 0, 0, 1, 1)
        self._axialControlsLayout.addWidget(self._axialButtons, 1, 0, 1, 1)
        
        self.addWidget(self._axialControls)
        
    def SetState(self, state):
        if "axialRoiPos" in state:
            self.tab.AxialClearRois()
            rois = state["axialRoiPos"]   
            for roiState in rois:
                roi = Utils.MyEllipseROI((0, 0), (1e-6, 1e-6))
                roi.setState(roiState)
                self.tab.AxialAddRoi(roi, samePos = True)
            
            if len(self.tab._axialRois) < 1:
                self.tab.AxialAddRoi(Utils.MyEllipseROI((0.0, 0.0), (2e-6, 2e-6), angle = 0.0))
            
    def SaveState(self):
        res = {}
        res["axialRoiPos"] = [roi.saveState() for roi in self.tab._axialRois]
        return res
        
    def DrawPlot(self):
        # profile = cProfile.Profile()
        # profile.enable()
        
        self.main._curPSF.UpdateCalibInterpolation()
        calibValues = self.main._curPSF.zCalibParamValues
        
        pens = ["b", "r", "g", "y"]
        symbols = ["x", "+", "o", "s"]
        for i, (name, calibParam) in enumerate(calibValues):
            if calibParam is None:
                continue
            
            for j, fitPixels in enumerate(calibParam.fitPixelsList):
                data = calibParam.GetData(fitPixels)
                zs, means, stds = data["zs"], data["means"], data["stds"]
                
                self._axialPlot.plot(zs, means, symbol = symbols[j], symbolPen = pens[i], \
                                     pen = None, name = "   %s (%d px)" % (name, fitPixels))
                
                if calibParam.useStds:
                    self._axialPlotErr = pg.ErrorBarItem(x = zs, y = means, top = stds, \
                                                         bottom = stds, pen = pens[i], beam = 4e-9)
                    self._axialPlot.plotItem.addItem(self._axialPlotErr)
                    
                if "interpCoefs" in data:
                    zsAxialInterp = np.linspace(zs[0], zs[-1], 500)
                    interpValues = calibParam.GetInterpolatedValue(zsAxialInterp, fitPixels)                    
                    self._axialPlot.plot(zsAxialInterp, interpValues, pen = pens[i])
                
        self._axialPlot.setLabels(left = ("", "m"), bottom = ("z", "m"))
        # profile.disable()
        # pstats.Stats(profile).sort_stats("cumulative").print_stats(15) 
    
    def ClearPlot(self):
        self._axialPlot.clear()
        self._axialPlot.plotItem.legend.items = []
    
    def Autoscale(self):
        pass
    
    def _CalibrateAxialPosition(self):
        self.main.ui.statusbar.showMessage("Axial calibration started...")
        self._calibButton.setEnabled(False)
        
        roisState = [roi.saveState() for roi in self.tab._axialRois]
        try:
            self.main._curAxialCalibrator.CalibratePsf(self.main._curPSF,
                                                  self.main._axialFseries, \
                                                  roisState)
            self.tab.DrawFrame()
        finally:
            self.main.ui.statusbar.showMessage("Axial calibration done.")
            self._calibButton.setEnabled(True)
            
    def _DeleteAxialCalibration(self):
        self.main._curPSF.SetAxialCalibPoints([])
        self.tab.DrawFrame()
        
#===============================================================================
# WobblePlotDock
#===============================================================================

class WobblePlotDock(Utils.PlotWidgetDock):
        
    def DrawPlot(self):
        self.main._curPSF.UpdateCalibInterpolation()
        calibValues = self.main._curPSF.wobbleCalibParamValues
        
        pens = ["b", "r", "g", "y"]
        symbols = ["x", "+", "o", "s"]
        for i, (name, calibParam) in enumerate(calibValues):
            if calibParam is None:
                continue
            
            for j, fitPixels in enumerate(calibParam.fitPixelsList):
                data = calibParam.GetData(fitPixels)
                zs, means, stds = data["zs"], data["means"], data["stds"]
                
                self.plt.plot(zs, means, symbol = symbols[j], symbolPen = pens[i], \
                                     pen = None, name = "   %s (%d px)" % (name, fitPixels))
                
                if calibParam.useStds:
                    self.pltErr = pg.ErrorBarItem(x = zs, y = means, top = stds, \
                                                         bottom = stds, pen = pens[i], beam = 4e-9)
                    self.plt.plotItem.addItem(self.pltErr)
                    
                if "interpCoefs" in data:
                    zsAxialInterp = np.linspace(zs[0], zs[-1], 500)
                    interpValues = calibParam.GetInterpolatedValue(zsAxialInterp, fitPixels)                    
                    self.plt.plot(zsAxialInterp, interpValues, pen = pens[i])
                
        self.plt.setLabels(left = ("", "m"), bottom = ("z", "m"))
  
#===============================================================================
# AxialCalibTab
#===============================================================================

class AxialCalibTab(Utils.DockAreaTabWidgetBase, Utils.MovableFramesBase):
    
    def __init__(self, *args, **kwargs):
        # Init super-classes
        Utils.DockAreaTabWidgetBase.__init__(self, *args, **kwargs)
        Utils.MovableFramesBase.__init__(self)
        self._InitDocks()
        
        # Init variables
        self._axialRois = []
        
    def _InitDocks(self):
        # Define docks        
        self._controlDock = CalibrationDock("Calibration", tab = self, main = self.main)
        self._wobbleDock = WobblePlotDock("Wobble", tab = self, main = self.main)
        self._originalDock = OriginalImageDock(tab = self, main = self.main)
        self._preprocessedDock = PreprocessedImageDock(tab = self, main = self.main)
        
        # Define default locations
        self._defaultDockPos = OrderedDict([
            (self._originalDock, ("left",)),
            (self._wobbleDock, ("right",)),
            (self._controlDock, ("above", self._wobbleDock)),
            (self._preprocessedDock, ("bottom", self._originalDock))])
        
        Utils.DockAreaTabWidgetBase._InitDocks(self)
        
    def GetImageCount(self):
        axialFseries = self.main._axialFseries
        if axialFseries is not None:
            return len(axialFseries)
        else:
            return 0
        
    def AxialAddRoi(self, oldRoi, samePos = False):
        pos = oldRoi.pos() if samePos else (0, 0)
        roi = Utils.MyEllipseROI(pos, oldRoi.size(), angle = oldRoi.angle())
        roi.sigRegionChangeFinished.connect(self._AxialRoiChangeFinished)
        roi.sigAddRequested.connect(self.AxialAddRoi)
        roi.sigRemoveRequested.connect(self.AxialRemoveRoi)
    
        self._originalDock.imv.addItem(roi)
        self._axialRois += [roi]
        
    def AxialRemoveRoi(self, roi, force = False):
        if len(self._axialRois) > 1 or force:
            self._axialRois.remove(roi)
            self._originalDock.imv.view.removeItem(roi)
        else:
            print "Cannot remove last ROI"
    
    def AxialClearRois(self):
        rois = self._axialRois[:]
        for roi in rois:
            self.AxialRemoveRoi(roi, force = True)
    
    def _AxialRoiChangeFinished(self):
        self.DrawFrame()


if __name__ == "__main__":
    pass
