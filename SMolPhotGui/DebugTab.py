import numpy as np
import scipy.stats
import pyqtgraph as pg
import pyqtgraph.exporters

from os import path
from PyQt4 import QtGui
from datetime import datetime
from collections import OrderedDict
from SMolPhot import Postprocessors, GroundTruthStats, RatingFuncGenerator

#===============================================================================
# Methods
#===============================================================================

PENS = ["b", "g", "r", "m", "y", "k", "c"]

def GetComparisonPlots(xlabel, plotNames = ["jac", "rmsXY", "rmsZ", "precision"]):
    res = OrderedDict()
    for plotName in plotNames:
        plt = pg.plot()
        plt.addLegend()
        plt.setTitle(str(plotName))
        plt.setLabels(bottom = (xlabel), left = (plotName))
        res[plotName] = plt
        
    if "rmsXY" in res:
        res["rmsXY"].setLabels(left = ("rmsXY", "m"))
    
    if "rmsZ" in res:
        res["rmsZ"].setLabels(left = ("rmsZ", "m"))
    return res

def SavePlots(plots, folder):
    for name, plot in plots.iteritems():
        exporter = pyqtgraph.exporters.ImageExporter(plot.plotItem)
        exporter.export(path.join(folder, "%s.png" % (name)))

#===============================================================================
# DebugTab
#===============================================================================

class DebugTab(QtGui.QScrollArea):
    
    def __init__(self, *args, **kwargs):
        self.main = kwargs.pop("main")
        tabName = kwargs.pop("tabName")
        QtGui.QScrollArea.__init__(self, *args, **kwargs)
        self.setFrameShape(QtGui.QFrame.NoFrame);
        
        self._widget = QtGui.QWidget()  
        self._layout = QtGui.QGridLayout(self._widget)
        
        self.setWidgetResizable(True)
        self.setWidget(self._widget)
        self.setObjectName(tabName)
        
        # Add controls
        
        self._buttonsWidget = QtGui.QWidget(self)
        self._buttonsLayout = QtGui.QVBoxLayout(self._buttonsWidget)
       
        self._compareGoodnessButton = QtGui.QPushButton("Compare goodness functions", self._buttonsWidget)
        self._compareGoodnessButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._compareGoodnessButton.pressed.connect(self._CompareGoodnessFunctions)
        
        self._compareErrorFuncsButton = QtGui.QPushButton("Compare error fucntions", self._buttonsWidget)
        self._compareErrorFuncsButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._compareErrorFuncsButton.pressed.connect(self._CompareErrorFunctions)
        
        self._compareBadnessFuncsButton = QtGui.QPushButton("Compare badness fucntions", self._buttonsWidget)
        self._compareBadnessFuncsButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._compareBadnessFuncsButton.pressed.connect(self._CompareBadnessFunctions)

        self._generateGoodnessFuncButton = QtGui.QPushButton("Generate goodness func", self._buttonsWidget)
        self._generateGoodnessFuncButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._generateGoodnessFuncButton.pressed.connect(self._GenerateGoodnessFunc)

        self._generateErrorFuncButton = QtGui.QPushButton("Generate error func", self._buttonsWidget)
        self._generateErrorFuncButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._generateErrorFuncButton.pressed.connect(self._GenerateErrorFunc)
        
        self._generateBadnessFuncButton = QtGui.QPushButton("Generate badness func", self._buttonsWidget)
        self._generateBadnessFuncButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._generateBadnessFuncButton.pressed.connect(self._GenerateBadnessFunc)
        
        self._compareCorrelationsButton = QtGui.QPushButton("Correlations", self._buttonsWidget)
        self._compareCorrelationsButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._compareCorrelationsButton.pressed.connect(self._CompareGoodnessCorrelation)
        
        self._gtStatsButton = QtGui.QPushButton("Gt stats", self._buttonsWidget)
        self._gtStatsButton.setIcon(QtGui.QIcon(":/16x16/16x16/Search_16x16.png"))
        self._gtStatsButton.pressed.connect(self._ShowGroundTruthStats)
        
        self._buttonsLayout.addWidget(self._compareGoodnessButton)
        self._buttonsLayout.addWidget(self._compareErrorFuncsButton)
        self._buttonsLayout.addWidget(self._compareBadnessFuncsButton)
        self._buttonsLayout.addWidget(self._generateGoodnessFuncButton)
        self._buttonsLayout.addWidget(self._generateErrorFuncButton)
        self._buttonsLayout.addWidget(self._generateBadnessFuncButton)
        self._buttonsLayout.addWidget(self._compareCorrelationsButton)
        self._buttonsLayout.addWidget(self._gtStatsButton)
        
        spacerItem = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Expanding, \
                                           QtGui.QSizePolicy.Minimum)
        self._buttonsLayout.addItem(spacerItem)
        
        self._layout.addWidget(self._buttonsWidget, 0, 0, 1, 1)
        spacerItem = QtGui.QSpacerItem(20, 20, QtGui.QSizePolicy.Expanding, \
                                       QtGui.QSizePolicy.Expanding)
        self._layout.addItem(spacerItem, 1, 0, 1, 1)
        
    def Shown(self):
        self._Update()
        
    def _Update(self):
        pass
    
    def _CompareGoodnessFunctions(self):
        goodnessFunctions = self.main._curPSF.goodnessFunctions.items()        
        patchPlots = GetComparisonPlots("# of patches", ["rmsXY", "rmsZ"])
        jacPlots = GetComparisonPlots("Jaccard (%)", ["rmsXY", "rmsZ"])
        locs = self.main._resultsTab._postProcessedLocs

        for i, (goodnessFuncName, goodnessFunc) in enumerate(goodnessFunctions):
            stat = GroundTruthStats(self.main._fseries, self.main._curGroundTruth)
            _, sortingRes = stat.SortedStatistics(locs, goodnessFunc)
            
            plotArgs = {"symbol": "x", "symbolPen": PENS[i], "pen": PENS[i], "name": goodnessFuncName}
            for paramName, plt in patchPlots.iteritems():
                plt.plot(sortingRes[paramName], **plotArgs)
            for paramName, plt in jacPlots.iteritems():
                plt.plot(sortingRes["jac"], sortingRes[paramName], **plotArgs)

    def _CompareErrorFunctions(self):
        errorFunctions = Postprocessors.ErrorRemovalPostprocessor.GetErrorFunctions()
        errorRemovalPP = self.main.GetPostprocessor(Postprocessors.ErrorRemovalPostprocessor)
        locs = errorRemovalPP._locsToProcess
    
        # Init rating func
        ratingFunc = RatingFuncGenerator(self.main._fseries, \
                                self.main._curGroundTruth, \
                                locs)
    
        plots = GetComparisonPlots("% to remove")
        # Loop over error functions
        for i, (name, errorFuncData) in enumerate(errorFunctions.iteritems()):
            r = ratingFunc.GetErrorFuncStats(errorFuncData)
            for paramName, plt in plots.iteritems():
                plt.plot(100.0 * r["partRemoved"], r[paramName], pen = PENS[i], \
                         symbol = "x", symbolPen = PENS[i], name = name)
        
    def _CompareBadnessFunctions(self):
        badnessFunctions = Postprocessors.BadLocsRemovalPostprocessor.GetBadnessFunctions()
        badLocsRemovalPP = self.main.GetPostprocessor(Postprocessors.BadLocsRemovalPostprocessor)
        locs = badLocsRemovalPP._locsToProcess
    
        # Init rating func
        ratingFunc = RatingFuncGenerator(self.main._fseries, \
                                self.main._curGroundTruth, \
                                locs)
    
        plots = GetComparisonPlots("% to remove")
        # Loop over error functions
        for i, (name, errorFuncData) in enumerate(badnessFunctions.iteritems()):
            r = ratingFunc.GetErrorFuncStats(errorFuncData)
            for paramName, plt in plots.iteritems():
                plt.plot(100.0 * r["partRemoved"], r[paramName], pen = PENS[i], \
                         symbol = "x", symbolPen = PENS[i], name = name)

    def _GenerateGoodnessFunc(self):
        saveParentFolder = "../../optimization/goodness"
        caseName = "%s %s" % (self.main._fseries.name, str(datetime.now()).replace(":", "-"))
        saveFolder = path.join(saveParentFolder, caseName)
        locs = self.main._resultsTab._resultsDock._allLocsNotPostprocessed
        
        rating = RatingFuncGenerator(self.main._fseries, self.main._curGroundTruth, locs)
        rating.GenerateGoodnessFunctions(saveFolder)
                
    
    def _GenerateErrorFunc(self):
        # Params
        partToRemove = 0.05
        nrOfIterationsList = [1, 2, 3, 5, 10, 20, 50]
        saveParentFolder = "../../optimization/error"
        caseName = "%s %s" % (self.main._fseries.name, str(datetime.now()).replace(":", "-"))
        saveFolder = path.join(saveParentFolder, caseName)
        bestnessCriteriaFunc = lambda x: (-x[2]["jac"], x[2]["rmsXY"])
        
        # Get locs
        errorRemovalPP = self.main.GetPostprocessor(Postprocessors.ErrorRemovalPostprocessor)
        locs = errorRemovalPP._locsToProcess
        
        # Generate func
        rating = RatingFuncGenerator(self.main._fseries, self.main._curGroundTruth, locs)
        info = rating.GenerateOptimalErrorFunction(partToRemove, nrOfIterationsList, \
                                                   bestnessCriteriaFunc, saveFolder)
 
        # Plotting
        plots = GetComparisonPlots("% to remove")
        for i, nrOfIterations in enumerate(nrOfIterationsList):
            iterationInfo = info[nrOfIterations]
            for paramName, plt in plots.iteritems():
                plt.plot(iterationInfo["partRemoved"], iterationInfo[paramName], \
                         symbol = "x", pen = PENS[i], symbolPen = PENS[i], \
                         name = "iter %d" % (nrOfIterations))
        SavePlots(plots, saveFolder)
        
    def _GenerateBadnessFunc(self):
        # Params
        partToRemove = 0.05
        nrOfIterationsList = [1, 2, 3, 5, 10, 20, 50]
        saveParentFolder = "../../optimization/badness"
        caseName = "%s %s" % (self.main._fseries.name, str(datetime.now()).replace(":", "-"))
        saveFolder = path.join(saveParentFolder, caseName)
        bestnessCriteriaFunc = lambda x: x[2]["rmsXY"]# ** 2.0 + x[2]["rmsZ"] ** 2.0
        
        # Get locs
        badLocsRemovalPP = self.main.GetPostprocessor(Postprocessors.BadLocsRemovalPostprocessor)
        locs = badLocsRemovalPP._locsToProcess
        
        # Generate func
        rating = RatingFuncGenerator(self.main._fseries, self.main._curGroundTruth, locs)
        info = rating.GenerateOptimalErrorFunction(partToRemove, nrOfIterationsList, \
                                                   bestnessCriteriaFunc, saveFolder)
        
        # Plotting
        plots = GetComparisonPlots("% to remove")
        for i, nrOfIterations in enumerate(nrOfIterationsList):
            iterationInfo = info[nrOfIterations]
            for paramName, plt in plots.iteritems():
                plt.plot(iterationInfo["partRemoved"], iterationInfo[paramName], \
                         symbol = "x", pen = PENS[i], symbolPen = PENS[i], \
                         name = "iter %d" % (nrOfIterations))
        SavePlots(plots, saveFolder)
            
    def _CompareGoodnessCorrelation(self):
        locs = self._postProcessedLocs[:]
        distXYs = np.array([loc.distXY for loc in locs])
        
        prop = {}
        prop["lsError"] = np.array([loc.lsError for loc in locs])
        prop["photons"] = np.array([loc.photons for loc in locs]).ravel()
        prop["amp"] = np.array([loc.amp for loc in locs]).ravel()
        prop["minFitDistXY"] = np.array([loc.minFitDistXY for loc in locs]).ravel()
        prop["sigmaX"] = np.array([loc.sigmaX for loc in locs]).ravel()
        prop["sigmaY"] = np.array([loc.sigmaY for loc in locs]).ravel()
        prop["offset"] = np.array([loc.offset for loc in locs]).ravel()
        prop["distXY"] = np.array([loc.distXY for loc in locs]).ravel()
        prop["fitStdAmp"] = np.array([loc.fitStdAmp for loc in locs]).ravel()
        prop["fitStdX0"] = np.array([loc.fitStdX0 for loc in locs]).ravel()
        prop["fitStdY0"] = np.array([loc.fitStdY0 for loc in locs]).ravel()
        prop["fitStdZ0"] = np.array([loc.fitStdZ0 for loc in locs]).ravel()
        prop["fitStdOffset"] = np.array([loc.fitStdOffset for loc in locs]).ravel()
        prop["pixels"] = np.array([float(loc.nrOfPixels) for loc in locs]).ravel()
    
        rValues = []
        for k, v in prop.iteritems():
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(distXYs, v)
            # print k, slope, intercept, "R2", r_value ** 2.0
            func = lambda d: slope * d + intercept
            rValues.append((k, slope, r_value ** 2.0))
            
            plt = pg.PlotItem()
            imv = pg.image(view = plt)
            
            plt.plot(distXYs, v, pen = None, symbol = "x", symbolSize = 2.0)
            plt.plot(distXYs, func(distXYs), pen = "r")
            plt.addLine(y = np.mean(v))
            plt.addLine(y = np.mean(v) + np.std(v))
            plt.addLine(y = np.mean(v) - np.std(v))
            
            plt.setLabels(bottom = ("distXY", "m"), left = (k))
            hist, xedges, yedges = np.histogram2d(distXYs, v, bins = 100)
            xStep, yStep = xedges[1] - xedges[0], yedges[1] - yedges[0]
            
            
            transform = QtGui.QTransform()
            transform.translate(xedges[0], yedges[0])
            transform.scale(xStep, yStep)
            imv.view.invertY(False)
            imv.view.setAspectLocked(False)    
            imv.setImage(hist, transform = transform)
            # plt.setRange(xRange = [0.0, 20e-9])
        
            
        print 
        print "Sorted R2"
        rValues.sort(key = lambda x:-x[-1])
        print "\n".join(["%s: slope %f, R2 %.3f" % (k, s, r2) for k, s, r2 in rValues])

    def _ShowGroundTruthStats(self):
        print "_ShowGroundTruthStats"    
        fseries = self.main._fseries
        # self._actNr = nrs[groundTruthToKeep].astype(int)
        # self._actFrameNr = frameNrs[groundTruthToKeep].astype(int) - 1
        # self._actCoords = 1e-9 * xs[groundTruthToKeep], 1e-9 * ys[groundTruthToKeep], 1e-9 * zs[groundTruthToKeep]
        # self._actIs = indensities[groundTruthToKeep]


        xs, ys, zs = fseries._actCoords
        
        pltXvsZ = pg.plot()
        pltXvsZ.setTitle("X vs Z")
        pltXvsZ.plot(xs, zs, symbol = "x", pen = None)
        pltXvsZ.setLabels(left = ("z", "m"), bottom = ("x", "m"))

        pltYvsZ = pg.plot()
        pltYvsZ.setTitle("Y vs Z")
        pltYvsZ.plot(ys, zs, symbol = "x", pen = None)
        pltYvsZ.setLabels(left = ("z", "m"), bottom = ("y", "m"))

        pltXvsY = pg.plot()
        pltXvsY.setTitle("X vs Y")
        pltXvsY.plot(xs, ys, symbol = "x", pen = None)
        pltXvsY.setLabels(left = ("y", "m"), bottom = ("x", "m"))
        
        # z cumulatiove distribution
        zCum = np.sort(zs)
        zProb = np.cumsum(np.ones_like(zCum)) / len(zCum)
        pltZProb = pg.plot()
        pltZProb.setTitle("Z vs prob")
        pltZProb.plot(zCum, 100 * zProb, symbol = "x", pen = None)
        pltZProb.setLabels(left = ("Cumulative probability", "%"), bottom = ("z", "m"))
        
        # intenisty cumulative distribution
        iCum = np.sort(fseries._actIs)
        iProb = np.cumsum(np.ones_like(iCum)) / len(iCum)
        pltIProb = pg.plot()
        pltIProb.setTitle("I vs prob")
        pltIProb.plot(iCum, 100 * iProb, symbol = "x", pen = None)
        pltIProb.setLabels(left = ("Cumulative probability", "%"), bottom = ("I"))

if __name__ == '__main__':
    pass
