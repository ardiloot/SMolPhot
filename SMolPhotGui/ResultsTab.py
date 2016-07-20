import Utils
import traceback
import numpy as np
import pyqtgraph as pg
import cProfile, pstats

from os import path
from time import time
from scipy import stats
from getpass import getuser
from PyQt4 import QtGui, QtCore
from collections import OrderedDict
from pyqtgraph.dockarea import Dock
from SMolPhot import GroundTruthStats, PublishToToplist, SaveResultsToFile
from SMolPhot.PSFs import DuplicateLocsList
     
#===============================================================================
# ResultDock
#===============================================================================

class ResultDock(Dock, Utils.SaveStateBaseClass):
    
    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = False
        Dock.__init__(self, *args, **kwargs)   
        self.plt = pg.PlotWidget()
        
        self._controlsWidget = QtGui.QWidget(self)
        self._controlsLayout = QtGui.QVBoxLayout(self._controlsWidget)
        self._buttonsWidget = QtGui.QWidget(self._controlsWidget)
        self._buttonsLayout = QtGui.QHBoxLayout(self._buttonsWidget)
        self._processFramesButton = QtGui.QPushButton("Process frames", self._buttonsWidget)
        self._processFramesButton.setIcon(QtGui.QIcon(":/16x16/16x16/Play_16x16.png"))
        self._processFramesButton.pressed.connect(self.tab._ProcessFramesPressed)
        self._framesFrom = QtGui.QSpinBox(self._buttonsWidget)
        self._framesFrom.setEnabled(False)
        self._framesTo = QtGui.QSpinBox(self._buttonsWidget)
        self._studyName = QtGui.QLineEdit(self._buttonsWidget)
        self._studyName.setText("test study")
        self._studyName.setFixedWidth(120)
        self._saveResultsButton = QtGui.QPushButton("Save", self._buttonsWidget)
        self._saveResultsButton.setIcon(QtGui.QIcon(":/16x16/16x16/Save_16x16.png"))
        self._saveResultsButton.pressed.connect(self.tab._SaveResults)
        
        self._buttonsLayout.addWidget(QtGui.QLabel("Study name: "))
        self._buttonsLayout.addWidget(self._studyName)
        self._buttonsLayout.addWidget(QtGui.QLabel("From: "))
        self._buttonsLayout.addWidget(self._framesFrom)
        self._buttonsLayout.addWidget(QtGui.QLabel("To: "))
        self._buttonsLayout.addWidget(self._framesTo)
        self._buttonsLayout.addWidget(self._processFramesButton)
        self._buttonsLayout.addWidget(self._saveResultsButton)
        spacerItem = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Expanding, \
                                           QtGui.QSizePolicy.Minimum)
        self._buttonsLayout.addItem(spacerItem)
        
        # Progress bar
        self._resultsProgressBar = QtGui.QProgressBar(self._controlsWidget)
        self._resultsProgressBar.setEnabled(True)
        self._resultsProgressBar.setProperty("value", 0)
        
        # Parent layout
        self._controlsLayout.addWidget(self._buttonsWidget)
        self._controlsLayout.addWidget(self.plt)
        self._controlsLayout.addWidget(self._resultsProgressBar)
        self.addWidget(self._controlsWidget)
        
        # Variables
        self._allLocsNotPostprocessed = []
            
    def SetState(self, state):
        if "studyName" in state:
            self._studyName.setText(state["studyName"])
        
    def SaveState(self):
        return {"studyName": str(self._studyName.text())}
    
    def DrawPlot(self):
        startTime = time()
        thread = self.tab._processFramesThread
        
        # Set process
        process = int(100.0 * thread.progress) 
        self._resultsProgressBar.setValue(process)

        # Plot
        newLocs = thread.GetNewLocs()
        self._allLocsNotPostprocessed += newLocs 
        xPoints = [l.x for l in newLocs]
        yPoints = [l.y for l in newLocs]
        
        self.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
        self.plt.plot(xPoints, yPoints, symbol = "x", pen = None, \
                           symbolPen = "r", symbolBrush = None, symbolSize = 2.0, \
                           antialias = False)
        
        timeTotal = time() - startTime
        print "_ResultsTabUpdatePlots time", timeTotal, len(xPoints), len(self._allLocsNotPostprocessed)
        
        # Stats
        if self.main._fseries.HasGroundTruth():
            statStr = thread.stats.GetStatsStr()
            self.main.ui.statusbar.showMessage("Runtime stats: %s" % (statStr))
            
        QtGui.QApplication.processEvents()
    
    def ClearPlot(self):
        self.plt.clear()
        
    def Autoscale(self):
        pass
        
#===============================================================================
# GoodnessDock
#===============================================================================

class GoodnessDock(Dock, Utils.SaveStateBaseClass):
    
    class PatchVsRmsXYDock(Utils.PlotWidgetDock):
    
        def DrawPlot(self):
            if self.tab.tab.tab._sortingResReal is None:
                return
            rmsXY = self.tab.tab.tab._sortingResReal["rmsXY"]
            #lastPatch = self.tab.tab.tab._sortingResReal["lastPatch"]
            pen = "b"
            self.plt.setLabels(left = ("rmsXY", "m"), bottom = ("# of patch"))
            self.plt.plot(rmsXY, symbol = "x", symbolPen = pen, pen = pen)
            #self.plt.addLine(x = lastPatch - 1)
            
    class PatchVsRmsZDock(Utils.PlotWidgetDock):
    
        def DrawPlot(self):
            if self.tab.tab.tab._sortingResReal is None:
                return
            rmsZ = self.tab.tab.tab._sortingResReal["rmsZ"]  
            #lastPatch = self.tab.tab.tab._sortingResReal["lastPatch"] 
            pen = "b"
            self.plt.setLabels(left = ("rmsZ", "m"), bottom = ("# of patch"))
            self.plt.plot(rmsZ, symbol = "x", symbolPen = pen, pen = pen)
            #self.plt.addLine(x = lastPatch - 1)
            
    class PatchVsJacDock(Utils.PlotWidgetDock):
    
        def DrawPlot(self):
            if self.tab.tab.tab._sortingResReal is None:
                return
            jac = self.tab.tab.tab._sortingResReal["jac"]
            #lastPatch = self.tab.tab.tab._sortingResReal["lastPatch"]
            pen = "b"
            self.plt.setLabels(left = ("Jaccard (%)", "m"), bottom = ("# of patch"))
            self.plt.plot(jac, symbol = "x", symbolPen = pen, pen = pen)
            #self.plt.addLine(x = lastPatch - 1)

    class PatchVsGoodness(Utils.PlotWidgetDock):
    
        def DrawPlot(self):
            if self.tab.tab.tab._sortingResReal is None:
                return
            
            goodness = np.squeeze([loc.goodnessValue for loc in self.tab.tab.tab._postProcessedLocs])
            # lastPatch = self.tab.tab.tab._sortingResReal["lastPatch"]
            pen = "b"
            self.plt.setLabels(bottom = ("#"), left = ("goodness"))
            self.plt.plot(goodness, symbol = "x", symbolPen = pen, pen = pen)
            # self.plt.addLine(x = lastPatch - 1)

    class JacVsRmsXYDock(Utils.PlotWidgetDock):
    
        def __init__(self, *args, **kwargs):
            Utils.PlotWidgetDock.__init__(self, *args, **kwargs)
            
            self._saveLineAction = QtGui.QAction("Save line", self.plt)
            self._clearLastLineAction = QtGui.QAction("Clear last line", self.plt)
            self._clearLineAction = QtGui.QAction("Clear lines", self.plt)
            self._saveLineAction.triggered.connect(self.SaveLine)
            self._clearLastLineAction.triggered.connect(self.ClearLastLine)
            self._clearLineAction.triggered.connect(self.ClearLine)
            
            menu = self.plt.plotItem.vb.menu
            menu.addAction(self._saveLineAction)
            menu.addAction(self._clearLastLineAction)
            menu.addAction(self._clearLineAction)
            
            self.lastData = None
            self.extraLinesToPlot = []
        
        def SaveLine(self):
            if self.lastData is None:
                return
            self.extraLinesToPlot.append(self.lastData)
            self.ClearPlot()
            self.DrawPlot()
            
        def ClearLastLine(self):
            if len(self.extraLinesToPlot) < 1:
                return
            del self.extraLinesToPlot[-1]
            self.ClearPlot()
            self.DrawPlot()
            
        def ClearLine(self):
            self.extraLinesToPlot = []
            self.ClearPlot()
            self.DrawPlot()

        def DrawPlot(self):
            if self.tab.tab.tab._sortingResReal is None:
                return
            
            # Current
            jac = self.tab.tab.tab._sortingResReal["jac"]
            rmsXY = self.tab.tab.tab._sortingResReal["rmsXY"]
            
            # Optimal
            jacOptimal = self.tab.tab.tab._sortingResOptimalXY["jac"]
            rmsXYOptimal = self.tab.tab.tab._sortingResOptimalXY["rmsXY"]
            
            #lastPatch = self.tab.tab.tab._sortingResReal["lastPatch"] 
            pen = "b"
            self.plt.setLabels(left = ("rmsXY", "m"), bottom = ("Jaccard (%)"))
            self.plt.plot(jac, rmsXY, symbol = "x", symbolPen = pen, \
                          pen = pen, name = "currnet")
            self.plt.plot(jacOptimal, rmsXYOptimal, symbol = "x", \
                          symbolPen = "y", pen = "y", name = "optimal")
            
            for i, toPlot in enumerate(self.extraLinesToPlot):
                self.plt.plot(*toPlot, symbol = "x", symbolPen = "r", \
                              pen = "r", name = "saved #%d" % (i + 1))
            self.lastData = (jac, rmsXY)
            
            #self.plt.addLine(x = jac[lastPatch - 1])
    
    class JacVsRmsZDock(JacVsRmsXYDock):
    
        def DrawPlot(self):
            if self.tab.tab.tab._sortingResReal is None:
                return
            
            # Current
            jac = self.tab.tab.tab._sortingResReal["jac"]
            rmsZ = self.tab.tab.tab._sortingResReal["rmsZ"]
            
            # Optimal
            jacOptimal = self.tab.tab.tab._sortingResOptimalZ["jac"]
            rmsZOptimal = self.tab.tab.tab._sortingResOptimalZ["rmsZ"]
            
            #lastPatch = self.tab.tab.tab._sortingResReal["lastPatch"] 
            pen = "b"
            self.plt.setLabels(left = ("rmsZ", "m"), bottom = ("Jaccard (%)"))
            self.plt.plot(jac, rmsZ, symbol = "x", symbolPen = pen, \
                          pen = pen, name = "current")
            self.plt.plot(jacOptimal, rmsZOptimal, symbol = "x", \
                          symbolPen = "y", pen = "y", name = "optimal")
            
            for i, toPlot in enumerate(self.extraLinesToPlot):
                self.plt.plot(*toPlot, symbol = "x", symbolPen = "r", \
                              pen = "r", name = "saved #%d" % (i + 1))
            self.lastData = (jac, rmsZ)
            
            #self.plt.addLine(x = jac[lastPatch - 1])
        
    class GoodnessTab(Utils.DockAreaTabWidgetBase):
        def __init__(self, *args, **kwargs):
            # Init super-classes
            self.tab = kwargs.pop("tab")
            Utils.DockAreaTabWidgetBase.__init__(self, *args, **kwargs)
            self._InitDocks()

        def _InitDocks(self):
            # Define docks        
            self._patchVsRmsXYDock = GoodnessDock.PatchVsRmsXYDock("patch nr vs rmsXY", tab = self, main = self.main)
            self._patchVsRmsZDock = GoodnessDock.PatchVsRmsZDock("patch nr vs rmsZ", tab = self, main = self.main)
            self._patchVsJacDock = GoodnessDock.PatchVsJacDock("patch nr vs jac", tab = self, main = self.main)
            self._patchVsGoodnessDock = GoodnessDock.PatchVsGoodness("patch nr vs goodness", tab = self, main = self.main)
            self._jacVsRmsXYDock = GoodnessDock.JacVsRmsXYDock("jac vs rmsXY", tab = self, main = self.main)
            self._jacVsRmsZDock = GoodnessDock.JacVsRmsZDock("jac vs rmsZ", tab = self, main = self.main)
        
            # Define default locations
            self._defaultDockPos = OrderedDict([
                (self._patchVsRmsXYDock, ("left",)),
                (self._jacVsRmsXYDock, ("right",)),
                (self._patchVsRmsZDock, ("bottom", self._patchVsRmsXYDock)),
                (self._patchVsJacDock, ("bottom", self._patchVsRmsZDock)),
                (self._patchVsGoodnessDock, ("bottom", self._patchVsJacDock)),
                (self._jacVsRmsZDock, ("bottom", self._jacVsRmsXYDock))])
            
            Utils.DockAreaTabWidgetBase._InitDocks(self)

    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = True
        Dock.__init__(self, *args, **kwargs) 
        
        self.dockArea = GoodnessDock.GoodnessTab(main = self.main, tab = self, tabName = "goodnessTab")  
        self.addWidget(self.dockArea)
        
    def DrawPlot(self, *args, **kwargs):
        self.dockArea.DrawFrame(*args, **kwargs)

    def ClearPlot(self, *args, **kwargs):
        self.dockArea.ClearPlots(*args, **kwargs)
        
    def Autoscale(self, *args, **kwargs):
        self.dockArea.AutoscalePlots(*args, **kwargs)
        
#===============================================================================
# PostprocessedImageDock
#===============================================================================

class PostprocessedImageDock(Dock, Utils.SaveStateBaseClass):
    
    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = False
        Dock.__init__(self, *args, **kwargs)   
        self.plt = pg.PlotWidget()
        
        self._controlsWidget = QtGui.QWidget(self)
        self._controlsLayout = QtGui.QGridLayout(self._controlsWidget)
        self._buttonsWidget = QtGui.QWidget(self._controlsWidget)
        self._buttonsLayout = QtGui.QHBoxLayout(self._buttonsWidget)
        self._processFramesButton = QtGui.QPushButton("Update", self._buttonsWidget)
        self._processFramesButton.setIcon(QtGui.QIcon(":/16x16/16x16/Play_16x16.png"))
        self._processFramesButton.pressed.connect(self._ProcessFramesPressed)
        self._buttonsLayout.addWidget(self._processFramesButton)
        spacerItem = QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Expanding, \
                                           QtGui.QSizePolicy.Minimum)
        self._buttonsLayout.addItem(spacerItem)
        self._controlsLayout.addWidget(self._buttonsWidget, 0, 0, 1, 1)
        self._controlsLayout.addWidget(self.plt, 1, 0, 1, 1)
        self.addWidget(self._controlsWidget)
            
    def _ProcessFramesPressed(self):
        self.DrawPlot()
            
    def DrawPlot(self):
        self.ClearPlot()
        
        locs = self.tab._postProcessedLocs
        
        xPoints = [l.x for l in locs]
        yPoints = [l.y for l in locs]
        
        self.plt.plot(xPoints, yPoints, symbol = "x", pen = None, \
                           symbolPen = "r", symbolBrush = None, symbolSize = 2.0, \
                           antialias = False)
        self.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
    
    def ClearPlot(self):
        self.plt.clear()
        
    def Autoscale(self):
        pass
             
          
#===============================================================================
# ResultsTab
#===============================================================================

class ResultsTab(Utils.DockAreaTabWidgetBase):
    
    def __init__(self, *args, **kwargs):
        # Init super-classes
        Utils.DockAreaTabWidgetBase.__init__(self, *args, **kwargs)
        self._InitDocks()
        
        # Variables
        self._lastSavePath = ""
        self._resHistory = [] 
        self._processFramesThread = _ProcessFramesThread(self.main)
        self._processFramesThread.started.connect(self._ProcessFramesThreadStarted)
        self._processFramesThread.finished.connect(self._ProcessFramesThreadFinished)
        
        self._resultsUpdateTimer = QtCore.QTimer()
        self._resultsUpdateTimer.timeout.connect(self._resultsDock.DrawPlot)
        
        self._sortingResReal = None
        self._sortingResWithoutPP = None
        self._sortingResOptimalXY = None
        self._sortingResOptimalZ = None
        self._postProcessedLocs = []

    def _InitDocks(self):
        # Define docks        
        self._resultsDock = ResultDock("Hi-res image", tab = self, main = self.main)
        self._goodnessDock = GoodnessDock("Goodness", tab = self, main = self.main)
        self._postprocessedImageDock = PostprocessedImageDock("Postprocessed image", tab = self, main = self.main)
        
        
        # Define default locations
        self._defaultDockPos = OrderedDict([
            (self._goodnessDock, ("left",)),
            (self._postprocessedImageDock, ("above", self._goodnessDock)),
            (self._resultsDock, ("above", self._postprocessedImageDock))
            ])
        Utils.DockAreaTabWidgetBase._InitDocks(self)

    def Shown(self):
        if self.main._fseries is not None:
            self._resultsDock._framesFrom.setRange(1, len(self.main._fseries))
            self._resultsDock._framesFrom.setValue(1)
            self._resultsDock._framesTo.setRange(1, len(self.main._fseries))
            self._resultsDock._framesTo.setValue(len(self.main._fseries))
        Utils.DockAreaTabWidgetBase.Shown(self)

    def _UpdateHiResImage(self):
        self._resultsDock.DrawPlot(clear = False)

    def _ProcessFramesThreadStarted(self):
        self._resultsDock._processFramesButton.setText("Stop")
        self._resultsDock._processFramesButton.setIcon(QtGui.QIcon(":/16x16/16x16/Stop_16x16.png"))
        self._resultsDock._processFramesButton.setEnabled(True)
        self._resultsDock._framesFrom.setEnabled(False)
        self._resultsDock._framesTo.setEnabled(False)
        self.main.ui.parametersTree.setEnabled(False)
        self._resultsDock._studyName.setEnabled(False)
        self._resultsDock._resultsProgressBar.setEnabled(True)
        self._resultsDock._resultsProgressBar.setValue(0)
        self._resultsDock.ClearPlot()
        self._resultsDock._allLocsNotPostprocessed = []
        self._sortingResReal = None
        self._postProcessedLocs = []
        self._processFramesStartTime = time()
        self.main.ui.statusbar.showMessage("Started processing frames...")
        self._resultsUpdateTimer.start(3000)
        
    def _ProcessFramesThreadFinished(self):
        self._resultsDock._processFramesButton.setText("Process frames")
        self._resultsDock._processFramesButton.setIcon(QtGui.QIcon(":/16x16/16x16/Play_16x16.png"))
        self._resultsDock._processFramesButton.setEnabled(True)
        self._resultsDock._framesFrom.setEnabled(False)
        self._resultsDock._framesTo.setEnabled(True)
        self.main.ui.parametersTree.setEnabled(True)
        self._resultsDock._studyName.setEnabled(True)
        self._resultsDock._resultsProgressBar.setEnabled(False)
        self._resultsUpdateTimer.stop()
        
        if self._processFramesThread.exception is not None:
            raise self._processFramesThread.exception
        
        self._resultsDock.DrawPlot()
        self._DoPostProcessing()
        self.DrawFrame()
        
        # Print status
        processTime = time() - self._processFramesStartTime
        print "Frames processed in %.2f s" % (processTime)
        
    def _DoPostProcessing(self):
        # profile = cProfile.Profile()
        print "_DoPostProcessing"
        if self._processFramesThread.progress < 0.9999:
            print "Stats: processing was stopped"
            return
            
        # PostProcessing
        # profile.enable()        
        postProcessedLocs = DuplicateLocsList(self._resultsDock._allLocsNotPostprocessed)
        for postprocessor in self.main._postprocessors:
            postProcessedLocs = postprocessor.Apply(self.main._fseries, self.main._curPSF, self.main._curGroundTruth, postProcessedLocs)
        self._postProcessedLocs = postProcessedLocs
        # profile.disable()
        
        # Stats
        if self.main._fseries.HasGroundTruth():
            originalStats = self._processFramesThread.stats
            postprocessedStats = GroundTruthStats(self.main._fseries, self.main._curGroundTruth, \
                                                  originalStats._nrGroundTruthMolecules)
            postprocessedStats.AddLocations(postProcessedLocs)
            
            # Results history
            thread = self._processFramesThread
            stats = postprocessedStats.GetStats()
            stats.update({"framesFrom": int(thread.frameFrom), \
                 "framesTo": int(thread.frameTo), \
                 "process": thread.progress})
            self._resHistory.append((stats["jac"], stats["recall"], stats["precision"], 1e9 * stats["rmsXY"], 1e9 * stats["rmsZ"], 1e9 * stats["averageDx"], 1e9 * stats["averageDy"], 1e9 * stats["averageDz"], stats["avgPhotonsScale"], stats["avgPhotonsDelta"]))
        
            print "Res history:"
            print "\t".join(["jac", "recall", "pres", "rmsXY", "rmsZ", "avgeDx", "avgDy", "avgDz", "phDelta", "phScale"])
            for hist in self._resHistory:
                print "\t".join(["%.3f" % (v) for v in hist])
                
            print "Without post-processing:"
            s = originalStats.GetStats()
            print "\t".join(["%.3f" % (v) for v in [s["jac"], s["recall"], s["precision"], 1e9 * s["rmsXY"], 1e9 * s["rmsZ"], 1e9 * s["averageDx"], 1e9 * s["averageDy"], 1e9 * s["averageDz"], s["avgPhotonsScale"], s["avgPhotonsDelta"]]])
            print
       
            # Real statistics
            statisticsReal = GroundTruthStats(self.main._fseries, self.main._curGroundTruth)
            _, self._sortingResReal = statisticsReal.SortedStatistics(postProcessedLocs, \
                lambda loc: loc.goodnessValue, alreadySorted = True)
       
            # Optimal statistics
            statOptimalXY = GroundTruthStats(self.main._fseries, self.main._curGroundTruth)
            _, self._sortingResOptimalXY = statOptimalXY.SortedStatistics(postProcessedLocs, lambda loc:-loc.distXY)
       
            statOptimalZ = GroundTruthStats(self.main._fseries, self.main._curGroundTruth)
            _, self._sortingResOptimalZ = statOptimalZ.SortedStatistics(postProcessedLocs, lambda loc:-loc.distZ)
    
            
            
            # Publish to toplist thread.frameFrom 
            toPublish = {"userName": getuser(),
                         "datasetName": self.main._fseries.name,
                         "testName": str(self._resultsDock._studyName.text()),
                         "framesFrom": thread.frameFrom,
                         "framesTo": thread.frameTo,
                         "goodnessFunc": "None",
                         "conf": self.main.GetConf(removeDockingStates = True),
                         "results": {"originalStats": originalStats.GetStats(),
                                     "postprocessedStats": postprocessedStats.GetStats(),
                                     "goodness": self._sortingResReal}}
        
            try:
                PublishToToplist(toPublish)
            except:
                traceback.print_exc()
        
            self.main.ui.statusbar.showMessage("Postprocessed stats: %s" % (postprocessedStats.GetStatsStr()))
        else:
            # No ground-truth
            self.main.ui.statusbar.showMessage("Processing done.")

        # print "Post-processing profile:"
        # pstats.Stats(profile).sort_stats("cumtime").print_stats(30)

    def _ProcessFramesPressed(self):
        self._resultsDock._processFramesButton.setEnabled(False)
        thread = self._processFramesThread
        
        if thread.isRunning():
            thread.Stop()
        else:
            thread.frameFrom = self._resultsDock._framesFrom.value() - 1
            thread.frameTo = self._resultsDock._framesTo.value()
            thread.start(priority = QtCore.QThread.HighPriority)
    
    def _SaveResults(self):        
        fname = str(QtGui.QFileDialog.getSaveFileName(self, "Save loacations", \
                                                       path.dirname(self._lastSavePath), \
                                                       "CSV (*.csv)"))
        
        if fname == "":
            return
        
        self._lastSavePath = fname
        SaveResultsToFile(self._postProcessedLocs, fname)
     
#===============================================================================
# _ProcessFramesThread
#===============================================================================
 
class _ProcessFramesThread(QtCore.QThread):
        
        def __init__(self, main, *args, **kwargs):
            QtCore.QThread.__init__(self, *args, **kwargs)
            self._mutex = QtCore.QMutex()
            self.main = main
            self.frameFrom, self.frameTo = 0, 0
            self._Init()
            
        def _Init(self):
            self.stopped = False
            self.progress = 0
            self.locs = []
            
            fseries = self.main._fseries
            if fseries is not None:
                gt = self.main._curGroundTruth
                if self.main._fseries.HasGroundTruth():
                    nrGroundTruthMolecules = fseries.GetNumberOfActMolecules(self.frameFrom, self.frameTo)
                    self.stats = GroundTruthStats(fseries, gt, nrGroundTruthMolecules)
                
        def run(self):
            self.exception = None
            try:
                profile = cProfile.Profile()
                psf = self.main._curPSF
                localizer = self.main._curLocalizer
                axialCalibrator = self.main._curAxialCalibrator
                fitMode = "z" if psf.hasCalibrationData else "sigma"
                hasGroundTruth = self.main._fseries.HasGroundTruth()
                
                for frameNr in range(self.frameFrom, self.frameTo):
                    if self.stopped:
                        return
                    
                    # Get frame
                    frame = self.main._fseries.GetPreprocessedFrame(frameNr)
                    
                    # Find molecules
                    profile.enable()
                    locsNew = localizer.FindMolecules(frame, psf, axialCalibrator, fitMode = fitMode)
                    profile.disable()
                    
                    # Add found locations
                    self._mutex.lock()
                    self.locs += locsNew
                    self._mutex.unlock()
                    
                    # Ground-truth
                    if hasGroundTruth:
                        self._mutex.lock()
                        self.stats.AddLocations(locsNew)
                        self._mutex.unlock()

                    # Update process
                    self.progress = float(frameNr - self.frameFrom) / (self.frameTo - self.frameFrom - 1)   
                
                pstats.Stats(profile).sort_stats("cumtime").print_stats(50)
            except Exception, ex:
                traceback.print_exc()
                self.exception = ex
                
        def start(self, *args, **kwargs):
            self._Init()
            return QtCore.QThread.start(self, *args, **kwargs)
  
        def Stop(self):
            self.stopped = True
            
        def GetNewLocs(self):
            self._mutex.lock()
            res = self.locs
            self.locs = []
            self._mutex.unlock()
            return res
   
    
if __name__ == '__main__':
    pass
