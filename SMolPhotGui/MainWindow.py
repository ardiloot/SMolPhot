import os
import yaml
import Utils
import appdirs
import SMolPhot
import traceback
import pyqtgraph as pg

from os import path
from datetime import datetime
from PyQt4 import QtCore, QtGui
from collections import OrderedDict
from QtDesigner import Ui_MainWindow
from SMolPhot import ModuleConf, FrameSeries
from pyqtgraph.parametertree import Parameter

from InfoTab import InfoTab
from AxialCalibTab import AxialCalibTab
from FramesTab import FramesBaseTab
from ResultsTab import ResultsTab
from DebugTab import DebugTab

APPNAME = "SMolPhot"
APPAUTHOR = "MolPhot"
DEFAULT_CONF_FILE = path.join(appdirs.user_config_dir(APPNAME, APPAUTHOR, roaming = True), \
                              "config.yaml")

#===============================================================================
# MainWindow
#===============================================================================
class MainWindow(QtGui.QMainWindow):
    
    def __init__(self, parent = None, developerMode = False):
        self._firstShow = True
        self._developerMode = developerMode
        super(MainWindow, self).__init__(parent)
        self.ui = Ui_MainWindow.Ui_MainWindow()
        self.ui.setupUi(self)
        
        # Window state
        self.settings = QtCore.QSettings(APPAUTHOR, APPNAME)
        self.restoreGeometry(self.settings.value("geometry").toByteArray())
        self.restoreState(self.settings.value("windowState").toByteArray())
        
        # Config pyQtGraph
        pg.setConfigOption("background", "w")
        pg.setConfigOption("foreground", "k")
        
        # Init backend
        self._fseries = None
        self._axialFseries = None
        self._InitModules()
        
        # Load configuration
        self._preloadFrameSeries = None
        self._saveToConf += [("GUI", self)]
        self._LoadConf()
        
        # Init GUI
        self._ConnectActions()
        self._InitTabs()
        self._InitParametersDock()
        
        # Disable
        self._GuiSetEnabled(False)
        
        
    # Overrides ----------------------------------------------------------------
        
    def showEvent(self, *args, **kwargs):
        if self._firstShow:
            self._firstShow = False
            self._CurrentTabChanged()
            if self._developerMode and self._preloadFrameSeries is not None:
                self._LoadFrameSeries(self._preloadFrameSeries, preload = True)
        return QtGui.QMainWindow.showEvent(self, *args, **kwargs)
    
    def closeEvent(self, event):
        self.settings.setValue("geometry", self.saveGeometry())
        self.settings.setValue("windowState", self.saveState())
        self._SaveConf()
        QtGui.QMainWindow.closeEvent(self, event)
    
    # Backend-------------------------------------------------------------------
        
    def _InitModules(self):
        (self._preprocessors, \
         self._localizers, \
         self._axialCalibrators, \
         self._psfs, \
         self._postprocessors, \
         self._groundTruths), self._saveToConf = ModuleConf.GetModuleConf()
         
        self._curPSF = self._psfs[0]
        self._curAxialCalibrator = self._axialCalibrators[0]
        self._curLocalizer = self._localizers[0]
        self._curGroundTruth = self._groundTruths[0]
        
    def _LoadFrameSeries(self, filename, preload = False):
        try:
            print "Loading frame series..."
            self._fseries = FrameSeries.FromMetafile(filename, \
                                                     self._preprocessors, \
                                                     series = "sequence",
                                                     maxFrameOverRide = None)
            print "Loading frame series done."
            
            print "Loading axial calibration..."
            self._axialFseries = FrameSeries.FromMetafile(filename, \
                                                          self._preprocessors, \
                                                          series = "axial calibration")
            print "Loading axial calibration done."
            
            # Optimize
            self._OptimizeForFrameSeries()
            
   
            
            # Axial calib tab        
            if preload and hasattr(self, "_preloadAxialCalibState"):
                self._axialCalibTab.SetState(self._preloadAxialCalibState)
            if len(self._axialCalibTab._axialRois) < 1:
                self._axialCalibTab.AxialAddRoi(Utils.MyEllipseROI((0.0, 0.0), \
                                                (2e-6, 2e-6), angle = 0.0))
    
            self._axialCalibTab.ChangeFrame(0)
            self._axialCalibTab.AutoscalePlots()
   
            # Frames tab
            if preload and hasattr(self, "_preloadFramesTabState"):
                self._framesTab.SetState(self._preloadFramesTabState)
            self._framesTab._framesTab.ChangeFrame(0)
            self._framesTab._iterationsTab.ChangeFrame(0)
            self._framesTab.AutoscalePlots()
            
   
            # Results tab
            if preload and hasattr(self, "_preloadResultsTabState"):
                self._resultsTab.SetState(self._preloadResultsTabState)
       
            self._CurrentTabChanged()
            self._GuiSetEnabled(True)
        except Exception, ex:
            ex.args += ("Error loading frame series %s" % (filename),)
            raise
        
    def _OptimizeForFrameSeries(self):
        self._fseries.OptimizePreprocessors()
        self._curLocalizer.OptimizeForSeries(self._fseries)
        
    def GetPostprocessor(self, postProcessorType):
        res = [pp for pp in self._postprocessors if type(pp) == postProcessorType]
        if len(res) <= 0:
            return None
        elif len(res) == 1:
            return res[0]
        else:
            return res
        
    # Config -------------------------------------------------------------------
    
    def SetParams(self, **kwargs):
    
        if "dataset metafile" in kwargs:
            if self._firstShow:
                self._preloadFrameSeries = kwargs["dataset metafile"]
    
        if "current localizer" in kwargs:
            index = Utils.GetIndex([t.name for t in self._localizers], kwargs["current localizer"])
            self._curLocalizer = self._localizers[index]
            
        if "current psf" in kwargs:
            index = Utils.GetIndex([t.name for t in self._psfs], kwargs["current psf"])
            self._curPSF = self._psfs[index]
            
        if "current axial calibrator" in kwargs:
            index = Utils.GetIndex([t.name for t in self._axialCalibrators], kwargs["current axial calibrator"])
            self._curAxialCalibrator = self._axialCalibrators[index]
        
        if "current ground truth" in kwargs:
            index = Utils.GetIndex([t.name for t in self._groundTruths], kwargs["current ground truth"])
            self._curGroundTruth = self._groundTruths[index]
            
        if "axial tab state" in kwargs:
            if self._firstShow:
                self._preloadAxialCalibState = kwargs["axial tab state"]
            else:
                self._axialCalibTab.SetState(kwargs["axial tab state"])
            
        if "frames tab state" in kwargs:
            if self._firstShow:
                self._preloadFramesTabState = kwargs["frames tab state"]
            else:
                self._framesTab.SetState(kwargs["frames tab state"])
            
        if "results tab state" in kwargs:
            if self._firstShow:
                self._preloadResultsTabState = kwargs["results tab state"]
            else:
                self._resultsTab.SetState(kwargs["results tab state"])
            
    def GetParams(self):
        metaDataFile = self._fseries.metadatafile if self._fseries is not None else None
        
        res = {"dataset metafile": metaDataFile, \
                "current localizer": self._curLocalizer.name,
                "current psf": self._curPSF.name,
                "current axial calibrator": self._curAxialCalibrator.name,
                "current ground truth": self._curGroundTruth.name,
                "axial tab state": self._axialCalibTab.SaveState(),
                "frames tab state": self._framesTab.SaveState(),
                "results tab state": self._resultsTab.SaveState(),
                "date": str(datetime.now())}
        return res
    
    def _LoadConf(self, filename = DEFAULT_CONF_FILE):
        conf = {}
        try:
            with file(filename, "r") as stream:
                conf = yaml.safe_load(stream)
        except IOError, ex:
            print "Error Loading conf", ex
            traceback.print_exc()
        finally:
            try:
                SMolPhot.Helpers.SetConf(self._saveToConf, conf)
            except:
                traceback.print_exc()
    
    def _SaveConf(self, filename = DEFAULT_CONF_FILE, removeDockingStates = False):
        conf = self.GetConf(removeDockingStates = removeDockingStates)
        
        if not path.isdir(path.dirname(filename)):
            os.makedirs(path.dirname(filename))
        
        with file(filename, "w") as stream:
            yaml.safe_dump(conf, stream, default_flow_style = False)
        self.ui.statusbar.showMessage("Saved.")
        
    def GetConf(self, removeDockingStates = False):
        conf = SMolPhot.Helpers.GetConf(self._saveToConf, removeDockingStates = removeDockingStates)
        return conf
        
    # GUI Actions --------------------------------------------------------------
        
    def _ConnectActions(self):
        self.ui.actionOpen.triggered.connect(self._OpenFrameSeries)
        self.ui.actionSave.triggered.connect(self._SaveButtonPressed)
        self.ui.actionExit.triggered.connect(self._Exit)
        
        self.ui.actionLoadConfFromFile.triggered.connect(self._LoadConfFromFile)
        self.ui.actionSaveConfToFile.triggered.connect(self._SaveConfToFile)
        
    def _OpenFrameSeries(self):
        directory = ""
        if self._fseries is None and self._preloadFrameSeries is not None:
            # First loading
            directory = path.dirname(self._preloadFrameSeries)
        
        fDialog = QtGui.QFileDialog(self, directory = directory)
        fDialog.setNameFilter("Datasets metafile (*.yml)");
        fDialog.setFileMode(QtGui.QFileDialog.ExistingFile)
        fDialog.setViewMode(QtGui.QFileDialog.Detail)
        if fDialog.exec_():
            filename = str(fDialog.selectedFiles()[0])
            self.ui.statusbar.showMessage("Loading file %s" % (filename))
            self._LoadFrameSeries(filename)
            self.ui.statusbar.showMessage("File %s loaded." % (filename))
            
    def _SaveButtonPressed(self):
        self._SaveConf()
                    
    def _LoadConfFromFile(self):
        fname = str(QtGui.QFileDialog.getOpenFileName(self, "Opens config", \
                                               "", \
                                               "YAML (*.yaml)"))
        if fname == "":
            return
        
        self._LoadConf(fname)
        self._InitParametersDock()
        
        self.ui.statusbar.showMessage("Config loaded from file.")
        
    def _SaveConfToFile(self):
        fname = str(QtGui.QFileDialog.getSaveFileName(self, "Save config", \
                                               "", \
                                               "YAML (*.yaml)"))
        if fname == "":
            return
        
        self._SaveConf(fname, removeDockingStates = True)    
        self.ui.statusbar.showMessage("Config saved to file.")
            
    def _Exit(self):
        self._SaveConf()
        self.close()
    
    
        
    # Tabs ---------------------------------------------------------------------
    
    def _InitTabs(self):
        self._tabWidget = QtGui.QTabWidget(self.ui.centralwidget)
        self._tabWidget.currentChanged.connect(self._CurrentTabChanged)
        self.ui.centralLayout.addWidget(self._tabWidget, 0, 0, 1, 1)
        
        self._infoTab = InfoTab(main = self, tabName = "infoTab")
        self._axialCalibTab = AxialCalibTab(main = self, tabName = "axialCalibTab")
        self._framesTab = FramesBaseTab(main = self, tabName = "framesTab")
        self._resultsTab = ResultsTab(main = self, tabName = "resultsTab")
        self._debugTab = DebugTab(main = self, tabName = "debugTab")
        
        self._tabWidget.addTab(self._infoTab, "Info")
        self._tabWidget.addTab(self._axialCalibTab, "Axial calibration")
        self._tabWidget.addTab(self._framesTab, "Frames")
        self._tabWidget.addTab(self._resultsTab, "Hi-res image")
        
        if self._developerMode:
            self._tabWidget.addTab(self._debugTab, "Debug")
        
    def _CurrentTabChanged(self):
        self.ui.statusbar.showMessage("")
        
        self._buttonsEnabled = False
        self._GetImageCount, self._DrawFrame, self._ClearPlots = \
            None, None, None
        self._AutoscalePlots, self._GetCurrentIndex, self._ChangeFrame = \
            None, None, None
        
        tabName = self._tabWidget.currentWidget().objectName()
        if tabName == "infoTab":
            self._infoTab.Shown()
        elif tabName == "axialCalibTab":
            self._axialCalibTab.Shown()
        elif tabName == "framesTab":
            self._framesTab.Shown()
        elif tabName == "resultsTab":
            self._resultsTab.Shown()
        elif tabName == "debugTab":
            self._debugTab.Shown()
        else:
            pass
    
    # Paramters tree -----------------------------------------------------------
    
    def _InitParametersDock(self):     
        
        # Preprocessors   
        self._paramsPreprocessors = Parameter.create(name = "Preprocessors", type = "group", \
            children = Utils.GetParamTreeDict([pp for pp in self._preprocessors], groupName = "Preprocessors", createGroup = False))
        
        # Localizers
        curLocalizer = self._localizers.index(self._curLocalizer)
        localizerDict = OrderedDict([(self._localizers[i].name, i) for i in range(len(self._localizers))])
        localizerSelector = [{"name": "Current localizer", "type": "list", "values": localizerDict, "value":  curLocalizer}]
        self._paramsLocalizer = Parameter.create(name = "Localizer", type = "group", children = localizerSelector)
        self._ChangeLocalizer(None)
        
        # PSFs
        curPsf = self._psfs.index(self._curPSF)
        psfsDict = OrderedDict([(self._psfs[i].name, i) for i in range(len(self._psfs))])
        psfSelector = [{"name": "Current PSF", "type": "list", "values": psfsDict, "value":  curPsf}]
        self._paramsPsf = Parameter.create(name = "PSF", type = "group", children = psfSelector)
        self._ChangePsf(None)

        # Axial calibrator
        curAxialCalibrator = self._axialCalibrators.index(self._curAxialCalibrator)
        axialCalibratorsDict = OrderedDict([(self._axialCalibrators[i].name, i) for i in range(len(self._axialCalibrators))])
        axialCalibratorsSelector = [{"name": "Current calibrator", "type": "list", "values": axialCalibratorsDict, "value":  curAxialCalibrator}]
        self._paramsAxialCalibrator = Parameter.create(name = "Axial calibrator", type = "group", children = axialCalibratorsSelector)
        self._ChangeAxialCalibrator(None)
        
        # Ground-truth
        curGroundTruth = self._groundTruths.index(self._curGroundTruth)
        groundTruthsDict = OrderedDict([(self._groundTruths[i].name, i) for i in range(len(self._groundTruths))])
        groundTruthSelector = [{"name": "Method", "type": "list", "values": groundTruthsDict, "value":  curGroundTruth}]
        self._paramsGroundTruth = Parameter.create(name = "Ground-Truth", type = "group", children = groundTruthSelector)
        self._ChangeGroundTruth(None)
        
        # Postprocessors
        self._paramsPostprocessors = Parameter.create(name = "Postprocessors", type = "group", \
            children = Utils.GetParamTreeDict([pp for pp in self._postprocessors], groupName = "Postprocessors", createGroup = False))
        
        self._parameters = [self._paramsPreprocessors,
                            self._paramsPsf,
                            self._paramsAxialCalibrator,
                            self._paramsLocalizer,
                            self._paramsGroundTruth,
                            self._paramsPostprocessors]
        
        # Add parameters
        self.ui.parametersTree.clear()
        for param in self._parameters:
            self.ui.parametersTree.addParameters(param, showTop = True)
            param.sigTreeStateChanged.connect(self._ParametersTreeChanged)
        self._paramsPostprocessors.sigTreeStateChanged.connect(self._ParametersTreeChangedPostProcessing)

    def _ParametersTreeChanged(self, param, changes):
        for param, change, data in changes:
            if change != "value":
                continue
            
            if "obj" in param.opts:
                obj, paramName = param.opts["obj"], param.opts["paramName"]
                obj.SetParams(**{paramName: data})
            elif param.opts["name"] == "Current localizer":
                self._ChangeLocalizer(data)
            elif param.opts["name"] == "Current PSF":
                self._ChangePsf(data)
            elif param.opts["name"] == "Method":
                self._ChangeGroundTruth(data)
            elif param.opts["name"] == "Current calibrator":
                self._ChangeAxialCalibrator(data)
            else:
                pass

            tabName = self._tabWidget.currentWidget().objectName()
            if tabName == "axialCalibTab":
                self._axialCalibTab.DrawFrame()
            elif tabName == "framesTab":
                self._framesTab.DrawFrame()
                
    def _ParametersTreeChangedPostProcessing(self, *args, **kwargs):
        tabName = self._tabWidget.currentWidget().objectName()
        if tabName == "resultsTab":
            self._resultsTab._DoPostProcessing()
        self._resultsTab.DrawFrame()
                
    def _ChangeLocalizer(self, nr):
        if nr != None:
            self._curLocalizer = self._localizers[nr]
            
        for ch in self._paramsLocalizer.childs[1:]:
            self._paramsLocalizer.removeChild(ch)
        self._paramsLocalizer.addChildren(Utils.GetParamTreeDict(self._curLocalizer))

    def _ChangePsf(self, nr):
        if nr != None:
            self._curPSF = self._psfs[nr]
            
        for ch in self._paramsPsf.childs[1:]:
            self._paramsPsf.removeChild(ch)
            
        self._paramsPsf.addChildren(Utils.GetParamTreeDict(self._curPSF))        

    def _ChangeAxialCalibrator(self, nr):
        if nr != None:
            self._curAxialCalibrator = self._axialCalibrators[nr]
            
        for ch in self._paramsAxialCalibrator.childs[1:]:
            self._paramsAxialCalibrator.removeChild(ch)
            
        self._paramsAxialCalibrator.addChildren(Utils.GetParamTreeDict(self._curAxialCalibrator))

    def _ChangeGroundTruth(self, nr):
        if nr != None:
            self._curGroundTruth = self._groundTruths[nr]
            
        for ch in self._paramsGroundTruth.childs[1:]:
            self._paramsGroundTruth.removeChild(ch)
        self._paramsGroundTruth.addChildren(Utils.GetParamTreeDict(self._curGroundTruth))     
        
    # GUI ----------------------------------------------------------------------
    
    def _GuiSetEnabled(self, enabled):
        self._tabWidget.setEnabled(enabled)
        self.ui.parametersTree.setEnabled(enabled)
        self.ui.parametersTree.setEnabled(enabled)

def Run():
    loggingDir = appdirs.user_log_dir(APPNAME, APPAUTHOR)
    SMolPhot.Helpers.SetupLogging(loggingDir = loggingDir)
    
    import sys
    app = QtGui.QApplication(sys.argv)
    main = MainWindow(developerMode = False)
    sys.excepthook = Utils.ExceptionHook
    main.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    Run()