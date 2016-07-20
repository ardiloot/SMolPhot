import Utils
import numpy as np
import cProfile, pstats

from PyQt4 import QtGui
from collections import OrderedDict
from SMolPhot import GroundTruthStats
from pyqtgraph.dockarea import Dock
from SMolPhot.PSFs import GetLocsCoordsArray

#===============================================================================
# FramesDock
#===============================================================================

class FramesDock(Dock, Utils.SaveStateBaseClass):

    class OriginalImageDock(Utils.ImageViewDock):
            
        def DrawPlot(self):
            fseries = self.main._fseries
            nr = self.tab._currentIndex
            
            # Load frame data
            originalFrame = fseries.GetOriginalFrame(nr)
    
            # Titles and labels
            self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
            self.label.setText("Orignial frame: #%d" % (nr + 1))
            
            # Actual plotting
            transform = QtGui.QTransform()
            transform.scale(fseries.pixelSize, fseries.pixelSize)    
            self.imv.setImage(originalFrame.data, autoLevels = False, \
                              autoRange = False, transform = transform)
    
        def name(self):
            return "Original"
        
    class PreprocessedImageDock(Utils.ImageViewDock):

        def DrawPlot(self):
            fseries = self.main._fseries
            nr = self.tab._currentIndex
            
            # Load frame data
            preprocessedFrame = self.tab._preprocessedFrame        
            
            # Titles and labels
            self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
            self.label.setText("Preprocessed frame: #%d" % (nr + 1))
         
            # Ground-truth        
            if preprocessedFrame.HasGroundTruth():
                gtXs, gtYs, _, _ = preprocessedFrame.GetActMolecules()
                self.imv.plt.plot(gtXs, gtYs, symbol = "+", pen = None, \
                                  symbolPen = "b", symbolBrush = None)
                
                idsTest, idsTrue, _, __ = self.tab._gtRes
                for idTest, idTrue in zip(idsTest, idsTrue):
                    testX, testY, _ = self.tab._locs[idTest].coord
                    trueX = fseries._actCoords[0][idTrue]
                    trueY = fseries._actCoords[1][idTrue]
                    # trueZ = fseries._actCoords[2][idTrue]
                    self.imv.plt.plot([testX, trueX], [testY, trueY], pen = "b")
         
            # Actual plotting
            transform = QtGui.QTransform()
            transform.scale(fseries.pixelSize, fseries.pixelSize)
            
            self.imv.plt.plot(*zip(*self.tab._molCoordsInitial), symbol = "o", \
                              pen = None, symbolPen = "r", symbolBrush = None)
            
            self.imv.plt.plot(*zip(*self.tab._molCoordsFit), symbol = "x", \
                              pen = None, symbolPen = "r", symbolBrush = None)
            
            self.imv.setImage(preprocessedFrame.data, autoLevels = False, \
                              autoRange = False, transform = transform)
            
        def name(self):
            return "Preprocessed"
    
    class SubtractedImageDock(Utils.ImageViewDock):
            
        def DrawPlot(self):
            fseries = self.main._fseries
            nr = self.tab._currentIndex
            
            # Titles and labels
            self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
            self.label.setText("Subtracted frame: #%d" % (nr + 1))
         
            if nr == 0:
                return
            
            # Load frame data
            preprocessedFrame = self.tab._preprocessedFrame
            lastPreprocessedFrame = fseries.GetPreprocessedFrame(nr - 1)        
            dif = preprocessedFrame.data - lastPreprocessedFrame.data

            # Actual plotting
            transform = QtGui.QTransform()
            transform.scale(fseries.pixelSize, fseries.pixelSize)
            self.imv.setImage(dif, autoLevels = False, \
                              autoRange = False, transform = transform)
            
        def name(self):
            return "Subtracted"

    class ResidualImageDock(Utils.ImageViewDock):
            
        def DrawPlot(self):
            # Load frame data
            fseries = self.main._fseries
            psf = self.main._curPSF
            nr = self.tab._currentIndex
            preprocessedFrame = self.tab._preprocessedFrame
                    
            # Residual
            residualFrame = preprocessedFrame.CopyFrame()
            for loc in self.tab._locs:
                psfData = psf.CalcWithoutOffset(preprocessedFrame.coords, *loc.bestFitParams)
                residualFrame.data -= psfData
            
            # Titles and labels
            self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
            self.label.setText("Residual: #%d" % (nr + 1))
         
            # Ground-truth        
            if preprocessedFrame.HasGroundTruth():
                gtXs, gtYs, _, _ = preprocessedFrame.GetActMolecules()
                self.imv.plt.plot(gtXs, gtYs, symbol = "+", pen = None, \
                                  symbolPen = "b", symbolBrush = None)
                
            # Actual plotting
            transform = QtGui.QTransform()
            transform.scale(fseries.pixelSize, fseries.pixelSize)
            self.imv.setImage(residualFrame.data, autoLevels = False, \
                              autoRange = False, transform = transform)
            
        def name(self):
            return "Residual"
    
    class BackgroundImageDock(Utils.ImageViewDock):
        
        def DrawPlot(self):
            # Load frame data
            fseries = self.main._fseries
            nr = self.tab._currentIndex
            preprocessedFrame = self.tab._preprocessedFrame
    
            # Titles and labels
            self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
            self.label.setText("Background: #%d" % (nr + 1))
         
            # Actual plotting
            transform = QtGui.QTransform()
            transform.scale(fseries.pixelSize, fseries.pixelSize)
            self.imv.setImage(preprocessedFrame._interpolatedBg, autoLevels = False, \
                              autoRange = False, transform = transform)
            
        def name(self):
            return "Background"
               
    class BgHistogramDock(Utils.PlotWidgetDock):
            
        def DrawPlot(self):
            preprocessedFrame = self.tab._preprocessedFrame
            
            # Hist
            hist, bin_edges = np.histogram(preprocessedFrame.data, 50, range = (-40.0, 40.0))
            binCenters = (bin_edges[1:] + bin_edges[:-1]) / 2.0
            self.plt.plot(binCenters, hist, symbol = "x", pen = None)
            
            # try:
            #    gauss = lambda x, A, mu, sigma: A * np.exp(-(x - mu) ** 2.0 / (2.0 * sigma ** 2.0))
            #    #fitVal, _ = optimize.curve_fit(gauss, binCenters, hist, [np.max(hist), 0.0, 100.0])
            #    binCentersFit = np.linspace(binCenters[0], binCenters[-1], 500)
            #    self.plt.plot(binCentersFit, gauss(binCentersFit, *fitVal))
            #    self.plt.addLine(x = fitVal[1])
            #    # self.plt.setLabels(left = ("counts"), bottom = ("z", "m"))
            #    print "Hist fit", fitVal  
            # except:
            #    pass
    
    class Tab(Utils.DockAreaTabWidgetBase, Utils.MovableFramesBase):
        def __init__(self, *args, **kwargs):
            # Init super-classes
            self.tab = kwargs.pop("tab")
            Utils.DockAreaTabWidgetBase.__init__(self, *args, **kwargs)
            Utils.MovableFramesBase.__init__(self)
            
            self._InitDocks()
            
        def _InitDocks(self):
            # Define docks
            self._originalDock = FramesDock.OriginalImageDock(tab = self, main = self.main)
            self._preprocessedDock = FramesDock.PreprocessedImageDock(tab = self, main = self.main)
            self._subtractedDock = FramesDock.SubtractedImageDock(tab = self, main = self.main)
            self._residualDock = FramesDock.ResidualImageDock(tab = self, main = self.main)
            # self._backgroundDock = FramesDock.BackgroundImageDock(tab = self, main = self.main)
            self._histDock = FramesDock.BgHistogramDock("Histogram", tab = self, main = self.main)
            
            
            # Define default locations
            self._defaultDockPos = OrderedDict([
                (self._originalDock, ("left",)),
                (self._preprocessedDock, ("right",)),
                (self._residualDock, ("bottom", self._originalDock)),
                # (self._backgroundDock, ("bottom", self._preprocessedDock)),
                (self._subtractedDock, ("above", self._residualDock)),
                (self._histDock, ("bottom", self._preprocessedDock))])
                        
            Utils.DockAreaTabWidgetBase._InitDocks(self)
            
        def GetImageCount(self):
            fseries = self.main._fseries
            if fseries is not None:
                return len(fseries)
            else:
                return 0
            
        def ChangeFrame(self, ind):
            Utils.MovableFramesBase.ChangeFrame(self, ind)
            # Update iterations tab
            self.tab.tab._iterationsTab.ChangeFrame(0)
            
            
        def DrawFrame(self, clear = True):
            profile = cProfile.Profile()
            profile.enable()
            
            # Init
            fseries, nr = self.main._fseries, self._currentIndex
            localizer = self.main._curLocalizer
            psf = self.main._curPSF
            axialCalibrator = self.main._curAxialCalibrator
            groundTruth = self.main._curGroundTruth
            preprocessedFrame = fseries.GetPreprocessedFrame(nr)  
            
            # Localize
            fitMode = "z" if psf.hasCalibrationData else "sigma"
            locs = localizer.FindMolecules(preprocessedFrame, psf, axialCalibrator, fitMode = fitMode)
            
            # Postprocessors
            #for postprocessor in self.main._postprocessors:
            #    locs = postprocessor.Apply(fseries, psf, groundTruth, locs)
    
            # Print debug info
            print "FRAME %d:" % (nr + 1)
            print "Found locations:"
            for loc in locs:
                print "%.3f %.3f %.3f, amp %.3f, offset %.3f" % \
                    (1e9 * loc.x, 1e9 * loc.y, 0.0 if loc.z is None else 1e9 * loc.z, loc.amp, loc.offset)
            print
                
            # Ground-truth        
            if preprocessedFrame.HasGroundTruth():
                # gtXs, gtYs, gtZs, gtIs = preprocessedFrame.GetActMolecules()
                # print "Ground-truth localtions:"
                # print "\n".join(["%.2f\t%.2f\t%.2f,\tI = %.2f" % \
                #    (1e9 * x, 1e9 * y, 1e9 * z, I) for x, y, z, I in zip(gtXs, gtYs, gtZs, gtIs)])
                
                gtStats = GroundTruthStats(self.main._fseries, self.main._curGroundTruth)
                self._gtRes = gtStats.MatchFrame(preprocessedFrame, locs)
                idsTest = self._gtRes[0]
                    
                # Failed locs
                for i in range(len(locs)):
                    if i in idsTest:
                        continue
                    print "Failed", 1e9 * locs[i].x, 1e9 * locs[i].y, \
                        (1e9 * locs[i].z if locs[i].z is not None else None), locs[i].bestFitParams
                
                # Stats
                if len(idsTest) > 0:
                    gtMolecules = len(fseries._actFrameMap[nr])
                    stats = GroundTruthStats(fseries, groundTruth, gtMolecules)
                    stats.AddLocations(locs)
                    print "Stats: ", stats.GetStatsStr(), "\n"
                    
            # Save for use in docks
            self._fitMode = fitMode
            self._locs = locs
            self._preprocessedFrame = preprocessedFrame
            self._molCoordsInitial = [(loc.ix, loc.iy) for loc in locs]
            self._molCoordsFit = [(loc.x, loc.y) for loc in locs]
            
            # Profile result
            print "Profile without drawing:"
            profile.disable()
            pstats.Stats(profile).sort_stats("cumtime").print_stats(20)
            
            # Draw
            Utils.DockAreaTabWidgetBase.DrawFrame(self, clear = clear)

    #--------------------------------------------------------------------------- 
    
    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = True
        
        Dock.__init__(self, *args, **kwargs)
        
        self._dockarea = FramesDock.Tab(main = self.main, tab = self, tabName = "FramesTab")  
        self.addWidget(self._dockarea)
        
    def DrawPlot(self, *args, **kwargs):
        self._dockarea.DrawFrame(*args, **kwargs)

    def ClearPlot(self, *args, **kwargs):
        self._dockarea.ClearPlots(*args, **kwargs)
        
    def Autoscale(self, *args, **kwargs):
        self._dockarea.AutoscalePlots(*args, **kwargs)
        
    def SetState(self, state):
        self._dockarea.SetState(state)

    def SaveState(self):
        return self._dockarea.SaveState()
        
#===============================================================================
# IterationsDock
#===============================================================================

class IterationsDock(Dock, Utils.SaveStateBaseClass):
    
    class StartFrameDock(Utils.ImageViewDock):
            
        def _DrawRegions(self, regions, **kwargs):
            fseries = self.main._fseries
            for region in regions:
                if region is None or region.bbox is None:
                    continue
                offset = np.array(region.offset)
                p0Index, p1Index = offset + region.bbox[:2], offset + region.bbox[2:]
                
                p0 = np.array(fseries.GetPixelCoords(p0Index)) - 0.5 * fseries.pixelSize
                p1 = np.array(fseries.GetPixelCoords(p1Index)) - 0.5 * fseries.pixelSize
                
                self.imv.plt.plot([p0[0], p1[0], p1[0], p0[0], p0[0]], \
                                   [p0[1], p0[1], p1[1], p1[1], p0[1]], **kwargs)
        
        def DrawPlot(self):
            if self.tab.GetImageCount() == 0:
                return
            
            # Init
            iterationNr = min(self.tab._currentIndex, self.tab.GetImageCount() - 1)
            iterationInfo = self.main._curLocalizer._iterationsInfoTmp[iterationNr]
            startFrame = iterationInfo["startFrame"]
            fseries = self.main._fseries
            
            # Titles and labels
            self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
            self.label.setText("Start frame: #%d / %d" % (iterationNr + 1, self.tab.GetImageCount()))
            
            # Actual plotting
            transform = QtGui.QTransform()
            transform.scale(fseries.pixelSize, fseries.pixelSize)    
            self.imv.setImage(startFrame.data, autoLevels = False, \
                              autoRange = False, transform = transform)
            
            # Locs
            if "locs" in iterationInfo:
                # Current iteration locs
                locsThis = iterationInfo["locs"]
                molCoordsInitial = [(loc.ix, loc.iy) for loc in locsThis]
                molCoordsFit = [(loc.x, loc.y) for loc in locsThis]
                self.imv.plt.plot(*zip(*molCoordsInitial), symbol = "o", pen = None, symbolPen = "r", symbolBrush = None)
                self.imv.plt.plot(*zip(*molCoordsFit), symbol = "x", pen = None, symbolPen = "r", symbolBrush = None)
        
                # Previous frame locs
                locPrevFrames = []
                for i in range(iterationNr):
                    locPrevFrames += self.main._curLocalizer._iterationsInfoTmp[i]["locs"]
                prevMolCoordsInitial = [(loc.ix, loc.iy) for loc in locPrevFrames]
                prevMolCoordsFit = [(loc.x, loc.y) for loc in locPrevFrames]
                self.imv.plt.plot(*zip(*prevMolCoordsInitial), symbol = "o", pen = None, symbolPen = "g", symbolBrush = None)
                self.imv.plt.plot(*zip(*prevMolCoordsFit), symbol = "x", pen = None, symbolPen = "g", symbolBrush = None)

                # Ground-truth        
                if fseries.HasGroundTruth():
                    locs = locsThis + locPrevFrames
                    gtXs, gtYs, _, __ = startFrame.GetActMolecules()
                    self.imv.plt.plot(gtXs, gtYs, symbol = "+", pen = None, symbolPen = "b", symbolBrush = None)
                    
                    gtStats = GroundTruthStats(self.main._fseries, self.main._curGroundTruth)
                    idsTest, idsTrue, _, __ = gtStats.MatchFrame(startFrame, locs)
                    for idTest, idTrue in zip(idsTest, idsTrue):
                        testX, testY, _ = locs[idTest].coord
                        trueX = fseries._actCoords[0][idTrue]
                        trueY = fseries._actCoords[1][idTrue]
                        #trueZ = fseries._actCoords[2][idTrue] 
                        self.imv.plt.plot([testX, trueX], [testY, trueY], pen = "b")
            
            # Regions   
            if "oldRegions" in iterationInfo:
                self._DrawRegions(iterationInfo["oldRegions"], symbol = None, pen = "g")
                    
            if "freshRegions" in iterationInfo:
                self._DrawRegions(iterationInfo["freshRegions"], symbol = None, pen = "y")
                
            if "updateRegions" in iterationInfo:
                self._DrawRegions(iterationInfo["updateRegions"], symbol = None, pen = "r")
                
            if "internalRegions" in iterationInfo:
                self._DrawRegions(iterationInfo["internalRegions"], symbol = None, pen = "g")
                
            if "freshRegions" in iterationInfo:
                self._DrawRegions(iterationInfo["freshRegions"], symbol = None, pen = "y")

            if "regions" in iterationInfo:
                self._DrawRegions(iterationInfo["regions"], symbol = None, pen = "r")
                                    
        def name(self):
            return "Start Frame"
        
    class PsfToDeleteDock(Utils.ImageViewDock):
                    
        def DrawPlot(self):
            if self.tab.GetImageCount() == 0:
                return
            
            # Init
            fseries = self.main._fseries
            iterationNr = self.tab._currentIndex
            iterationInfo = self.main._curLocalizer._iterationsInfoTmp[iterationNr]
            
            if "psfsToDelete" in iterationInfo:
                psfsToDelete = iterationInfo["psfsToDelete"]
                
                # Titles and labels
                self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
                self.label.setText("Psf to delete: #%d / %d" % (iterationNr + 1, self.tab.GetImageCount()))
                
                # Actual plotting
                transform = QtGui.QTransform()
                transform.scale(fseries.pixelSize, fseries.pixelSize)    
                self.imv.setImage(psfsToDelete, autoLevels = False, \
                                  autoRange = False, transform = transform)
            else:
                self.label.setText("Psf to delete: no data")
        
        def name(self):
            return "Psf to delete"
        
    class AfterDeleteDock(Utils.ImageViewDock):
            
        def DrawPlot(self):
            if self.tab.GetImageCount() == 0:
                return
            
            # Init
            fseries = self.main._fseries
            iterationNr = self.tab._currentIndex
            iterationInfo = self.main._curLocalizer._iterationsInfoTmp[iterationNr]

            if "frameNext" in iterationInfo:
                frameNext = iterationInfo["frameNext"]
                
                # Titles and labels
                self.imv.plt.setLabels(left = ("y", "m"), bottom = ("x", "m"))
                self.label.setText("After deletion: #%d / %d" % (iterationNr + 1, self.tab.GetImageCount()))
                
                # Actual plotting
                transform = QtGui.QTransform()
                transform.scale(fseries.pixelSize, fseries.pixelSize)    
                self.imv.setImage(frameNext.data, autoLevels = False, \
                                  autoRange = False, transform = transform)
            else:
                self.label.setText("After deletion: no data")
    
        def name(self):
            return "After delete"

    class Tab(Utils.DockAreaTabWidgetBase, Utils.MovableFramesBase):
        def __init__(self, *args, **kwargs):
            # Init super-classes
            self.tab = kwargs.pop("tab")
            Utils.DockAreaTabWidgetBase.__init__(self, *args, **kwargs)
            Utils.MovableFramesBase.__init__(self)
            
            self._InitDocks()
            
        def _InitDocks(self):
            # Define docks
            self._startFrameDock = IterationsDock.StartFrameDock(tab = self, main = self.main)
            self._psfsToDeleteDock = IterationsDock.PsfToDeleteDock(tab = self, main = self.main)
            self._afterDeleteDock = IterationsDock.AfterDeleteDock(tab = self, main = self.main)

            # Define default locations
            self._defaultDockPos = OrderedDict([
                (self._startFrameDock, ("left",)),
                (self._psfsToDeleteDock, ("right",)),
                (self._afterDeleteDock, ("bottom", self._startFrameDock))
                ])
            
            Utils.DockAreaTabWidgetBase._InitDocks(self)
            
        def GetImageCount(self):
            curLocalizer = self.main._curLocalizer
            if hasattr(curLocalizer, "_iterationsInfoTmp"):
                return len(curLocalizer._iterationsInfoTmp)
            else:
                return 0
            
        def DrawFrame(self, clear = True):
            Utils.DockAreaTabWidgetBase.DrawFrame(self, clear = clear)
            self.AutoscalePlots()

    #--------------------------------------------------------------------------- 

    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = True
        
        Dock.__init__(self, *args, **kwargs)
        
        self._dockarea = IterationsDock.Tab(main = self.main, tab = self, tabName = "IterationsTab")  
        self.addWidget(self._dockarea)
        
    def DrawPlot(self, *args, **kwargs):
        self._dockarea.DrawFrame(*args, **kwargs)

    def ClearPlot(self, *args, **kwargs):
        self._dockarea.ClearPlots(*args, **kwargs)
        
    def Autoscale(self, *args, **kwargs):
        self._dockarea.AutoscalePlots(*args, **kwargs)
        
    def SetState(self, state):
        self._dockarea.SetState(state)

    def SaveState(self):
        return self._dockarea.SaveState()
             
#===============================================================================
# FramesBaseTab
#===============================================================================

class FramesBaseTab(Utils.DockAreaTabWidgetBase):
    
    def __init__(self, *args, **kwargs):
        # Init super-classes
        Utils.DockAreaTabWidgetBase.__init__(self, *args, **kwargs)
        
        # Define docks
        self._framesDock = FramesDock("Frames", tab = self, main = self.main) 
        self._iterationsDock = IterationsDock("Iterations", tab = self, main = self.main)
        
        self._framesTab = self._framesDock._dockarea
        self._iterationsTab = self._iterationsDock._dockarea
        
        # Define default locations
        self._defaultDockPos = OrderedDict([
            (self._iterationsDock, ("left",)),
            (self._framesDock, ("above", self._iterationsDock))
            ])
        
        self._InitDocks()
        self._plotDocks.reverse()

if __name__ == '__main__':
    pass
