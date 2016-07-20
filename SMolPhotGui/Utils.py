import logging
import traceback
import numpy as np
import pyqtgraph as pg

from PyQt4 import QtGui, QtCore
from pyqtgraph import parametertree
from pyqtgraph.dockarea import DockArea, Dock

__all__ = []

logger = logging.getLogger(__name__)

#===============================================================================
# SaveStateBaseClass
#===============================================================================

class SaveStateBaseClass(object):
            
    def SetState(self, state):
        if len(state) > 0:
            raise NotImplementedError()
    
    def SaveState(self):
        return {}
    
#===============================================================================
# MovableFramesBase
#===============================================================================

class MovableFramesBase():
    
    def __init__(self, *args, **kwargs):
        # Variables for frame selection
        # Code for frame selection mainly copied form PyQtGraph code
        self._currentIndex = 0
        self._buttonsEnabled = True
        self._playTimer = QtCore.QTimer()
        self._playTimer.timeout.connect(self._PlayTimerTimeout)
        self._lastPlayTime = 0.0
        self._playRate = 0
        self._keysPressed = {}
        self._noRepeatKeys = [QtCore.Qt.Key_Left, QtCore.Qt.Key_Right, \
                              QtCore.Qt.Key_Up, QtCore.Qt.Key_Down, \
                              QtCore.Qt.Key_PageUp, QtCore.Qt.Key_PageDown]
        
    def keyPressEvent(self, ev):
        if self._buttonsEnabled and ev.key() in self._noRepeatKeys:
            ev.accept()
            if ev.isAutoRepeat():
                return
            self._keysPressed[ev.key()] = 1
            self._EvalKeyState()
        else:
            QtGui.QWidget.keyPressEvent(self, ev)
    
    def keyReleaseEvent(self, ev):
        if self._buttonsEnabled and ev.key() in self._noRepeatKeys:
            ev.accept()
            if ev.isAutoRepeat():
                return
            try:
                del self._keysPressed[ev.key()]
            except:
                self.keysPressed = {}
            self._EvalKeyState()
        else:
            QtGui.QWidget.keyReleaseEvent(self, ev)
            
    def ChangeFrame(self, ind):
        self._currentIndex = ind
        self.DrawFrame()
        
    def GetCurrentIndex(self):
        return self._currentIndex

    def GetImageCount(self):
        raise NotImplementedError()
        
    def DrawFrame(self, clear = True):
        raise NotImplementedError()

    def _EvalKeyState(self):        
        if len(self._keysPressed) == 1:
            key = list(self._keysPressed.keys())[0]
            if key == QtCore.Qt.Key_Right:
                self._Play(20)
                self._JumpFrames(1)
                self._lastPlayTime = pg.ptime.time() + 0.2
            elif key == QtCore.Qt.Key_Left:
                self._Play(-20)
                self._JumpFrames(-1)
                self._lastPlayTime = pg.ptime.time() + 0.2
            elif key == QtCore.Qt.Key_Up:
                self._Play(-100)
            elif key == QtCore.Qt.Key_Down:
                self._Play(100)
            elif key == QtCore.Qt.Key_PageUp:
                self._Play(-1000)
            elif key == QtCore.Qt.Key_PageDown:
                self._Play(1000)
        else:
            self._Play(0)
        
    def _Play(self, rate):
        self._playRate = rate
        if rate == 0:
            self._playTimer.stop()
            return
        
        self._lastPlayTime = pg.ptime.time()
        if not self._playTimer.isActive():
            self._playTimer.start(16)

    def _JumpFrames(self, n):
        ind = np.clip(self.GetCurrentIndex() + n, 0, self.GetImageCount() - 1)
        
        if ind == self.GetCurrentIndex():
            return
        
        self.ChangeFrame(ind)
               
    def _PlayTimerTimeout(self):
        now = pg.ptime.time()
        dt = now - self._lastPlayTime
        if dt < 0:
            return
        
        n = int(self._playRate * dt)
        if n != 0:
            self._lastPlayTime += (float(n) / self._playRate)
            if self.GetCurrentIndex() + n > self.GetImageCount():
                self._Play(0)
            self._JumpFrames(n)
        else:
            pass
                
        QtGui.QApplication.processEvents()


#===============================================================================
# MyImageView
#===============================================================================

class MyImageView(pg.ImageView):
    def __init__(self, *args, **kwargs):
        self.plt = pg.PlotItem()
        
        pg.ImageView.__init__(self, view = self.plt, *args, **kwargs)
        self.ui.menuBtn.setVisible(False)
        self.ui.roiBtn.setVisible(False)
        self.plt.invertY(False)
        
    def keyPressEvent(self, ev):
        QtGui.QWidget.keyPressEvent(self, ev)
    
    def keyReleaseEvent(self, ev):
        QtGui.QWidget.keyReleaseEvent(self, ev)
        
    def Clear(self):
        self.clear()
        self.plt.clearPlots()

#===============================================================================
# MyEllipseROI
#===============================================================================

class MyEllipseROI(pg.EllipseROI):
    sigAddRequested = QtCore.Signal(object)  # @UndefinedVariable
    
    def __init__(self, pos, size, **args):
        args["removable"] = True
        pg.EllipseROI.__init__(self, pos, size, **args)
        
        
    def paint(self, p, opt, widget):
        pg.EllipseROI.paint(self, p, opt, widget)
        
        p.setPen(pg.mkPen("r"))
        p.drawLine(QtCore.QLineF(0.5, 0.25, 0.5, 0.75))
        p.drawLine(QtCore.QLineF(0.25, 0.5, 0.75, 0.5))
        
    def getMenu(self):
        if self.menu is None:
            self.menu = QtGui.QMenu()
            self.menu.setTitle("ROI")
            
            addAct = QtGui.QAction("Add ROI", self.menu)
            addAct.triggered.connect(self.addRoiClicked)
            self.menu.addAction(addAct)
            self.menu.addAct = addAct
           
            remAct = QtGui.QAction("Remove ROI", self.menu)
            remAct.triggered.connect(self.removeClicked)
            self.menu.addAction(remAct)
            self.menu.remAct = remAct 
        return self.menu

    def addRoiClicked(self):
        # # Send remove event only after we have exited the menu event handler
        QtCore.QTimer.singleShot(0, lambda: self.sigAddRequested.emit(self))
        
    def saveState(self):
        state = pg.EllipseROI.saveState(self)
        state["angle"] = float(state["angle"])
        return state
     
#===============================================================================
# DockAreaTabWidgetBase
#===============================================================================

class DockAreaTabWidgetBase(QtGui.QWidget):
    
    def __init__(self, *args, **kwargs):
        self.main = kwargs.pop("main")
        tabName = kwargs.pop("tabName")
        QtGui.QWidget.__init__(self, *args, **kwargs)
        self.setObjectName(tabName)
        self._layout = QtGui.QGridLayout(self)
            
    def _InitDocks(self):
        
        # Define docking area
        if hasattr(self, "_dockArea"):
            self._dockArea.setParent(None)
        self._dockArea = DockArea()

        self._plotDocks = self._defaultDockPos.keys()
        
        # Add dock to area
        for dock, pos in self._defaultDockPos.iteritems():
            self._dockArea.addDock(dock, *pos)
        
        self._layout.addWidget(self._dockArea, 0, 0, 1, 1)
    
    def Shown(self):
        self.DrawFrame()
    
    def DrawFrame(self, clear = True):   
        if clear:
            self.ClearPlots()
            
        for dock in self._plotDocks:
            if not dock.automaticDraw:
                continue
            dock.DrawPlot()
        
    def ClearPlots(self):
        for dock in self._plotDocks:
            if not dock.automaticDraw:
                continue
            dock.ClearPlot()
            
    def AutoscalePlots(self):
        for dock in self._plotDocks:
            if not dock.automaticDraw:
                continue
            dock.Autoscale()
            
    def SaveState(self):
        res = {}
        res["dockingState"] = self._dockArea.saveState()
        for dock in self._plotDocks:
            res["dock_" + dock.name()] = dock.SaveState()
        return res
    
    def SetState(self, state):
        try:    
            if "dockingState" in state:
                self._dockArea.restoreState(state["dockingState"])           
        except:
            print "Docking area restore failed, restoring defaults:"
            traceback.print_exc()
            print "Restore defaults"
            self._InitDocks()
 
        for dock in self._plotDocks:
            stateName = "dock_" + dock.name()
            if stateName in state:
                dock.SetState(state[stateName])
            
 
#===============================================================================
# ImageViewDock
#===============================================================================

class ImageViewDock(Dock, SaveStateBaseClass):
    
    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = True
        
        Dock.__init__(self, self.name(), **kwargs)   
        self.imv = MyImageView()
        self.imv.setContentsMargins(2, 2, 2, 2)
        self.addWidget(self.imv)
        
    def DrawPlot(self):
        raise NotImplementedError()
    
    def ClearPlot(self):
        self.imv.Clear()
        
    def Autoscale(self):
        self.imv.autoLevels()

#===============================================================================
# PlotWidgetDock
#===============================================================================

class PlotWidgetDock(Dock, SaveStateBaseClass):
    
    def __init__(self, *args, **kwargs):
        self.tab = kwargs.pop("tab", None)
        self.main = kwargs.pop("main", None)
        self.automaticDraw = True
        
        Dock.__init__(self, *args, **kwargs)   
        self.plt = pg.PlotWidget()
        self.plt.addLegend()
        
        wrapWidget = QtGui.QWidget(self)
        wrapLayout = QtGui.QGridLayout(wrapWidget)
        wrapLayout.addWidget(self.plt, 0, 0, 1, 1)
        wrapWidget.setContentsMargins(2, 2, 2, 2)
        self.addWidget(wrapWidget)
        
    def DrawPlot(self):
        pass
    
    def ClearPlot(self):
        self.plt.clear()
        self.plt.plotItem.legend.items = []
        
    def Autoscale(self):
        pass

#===============================================================================
# CheckBoxSpinBox
#===============================================================================

class CheckBoxSpinBox(QtGui.QWidget):
    sigChanged = QtCore.Signal(object)
    sigChanging = QtCore.Signal(object, object)
    
    def __init__(self, *args, **kwargs):
        QtGui.QWidget.__init__(self, *args, **kwargs)
        
        
        self.checkbox = QtGui.QCheckBox()
        self.checkbox.toggled.connect(self.Changed)
        
        self.spinbox = pg.SpinBox()
        self.spinbox.sigValueChanged.connect(self.Changed)
        self.spinbox.sigValueChanging.connect(self.Changing)

        self.layout = QtGui.QHBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.layout.addWidget(self.checkbox)
        self.layout.addWidget(self.spinbox)
        
    def setValue(self, val):
        self.checkbox.setChecked(val["checked"])
        self.spinbox.setValue(val["value"])
        
    def value(self):
        return {"checked": self.checkbox.isChecked(), \
                "value": self.spinbox.value()}
        
    def setSpinBoxOpts(self, **opts):
        self.spinbox.setOpts(**opts)
        
    def Changed(self, obj):
        self.sigChanged.emit(obj)
        
    def Changing(self, obj, value):
        self.sigChanging.emit(obj, value)

#===============================================================================
# MyWidgetParameterItem
#===============================================================================

class MyWidgetParameterItem(parametertree.parameterTypes.WidgetParameterItem):

    def makeWidget(self):
        opts = self.param.opts
        t = opts['type']
        if t == '(bool, float)':
            defsSpinBox = {
                'min': None, 'max': None, 
                'step': 1.0, 'dec': False, 
                'siPrefix': False, 'suffix': ''
            }
            
            for k, v in opts.iteritems():
                if k in defsSpinBox:
                    defsSpinBox[k] = v
            defsSpinBox["value"] = opts["value"]["value"]
            
            if 'limits' in opts:
                defsSpinBox['bounds'] = opts['limits']
                
            w = CheckBoxSpinBox()
            w.setSpinBoxOpts(**defsSpinBox)
            self.hideWidget = False
        else:
            raise Exception("Unknown type '%s'" % unicode(t))
        return w
        
#===============================================================================
# TupleParameter
#===============================================================================

class TupleParameter(parametertree.Parameter):
    itemClass = MyWidgetParameterItem
    
    def __init__(self, *args, **kargs):
        parametertree.Parameter.__init__(self, *args, **kargs)
        if self.opts['type'] == 'color':
            self.value = self.colorValue
    
    def colorValue(self):
        return pg.mkColor(parametertree.Parameter.value(self))

parametertree.registerParameterType("(bool, float)", TupleParameter)

#===============================================================================
# Methods
#===============================================================================

def ExceptionHook(eType, value, tback):
    logger.error("ExceptionHook error!", exc_info = (eType, value, tback))
    
    box = QtGui.QMessageBox()
    box.setWindowTitle("Error")
    box.setIcon(QtGui.QMessageBox.Critical)
    
    box.setText("<b>Errors:</b> <br>  - " + "<br>   - ".join(reversed(map(str, value.args))))
    box.setInformativeText("Please contract your system administrator and show this message. \n" + \
                           "The application can not get started without a correct configuration.")
    box.setDetailedText("".join(traceback.format_exception(type, value, tback)))
    horizontalSpacer = QtGui.QSpacerItem(400, 0, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
    layout = box.layout()
    layout.addItem(horizontalSpacer, layout.rowCount(), 0, 1, layout.columnCount());
    box.exec_()
    
    # then call the default handler 
    #sys.__excepthook__(type, value, tback)  

def GetParamType(value):
    if type(value) == float or type(value) == np.float64:
        return "float"
    elif type(value) == int:
        return "int"
    elif type(value) == bool:
        return "bool"
    elif type(value) == str:
        return "str"
    elif type(value) == list:
        return "list"
    elif type(value) == dict:
        if "checked" in value and "value" in value:
            return "(bool, float)"
        else:
            raise ValueError("Unknown type %s, %s" % (type(value), value))
    else:
        raise ValueError("Unknown type %s, %s" % (type(value), value))

def GetParamTreeDict(obj, groupName = None, createGroup = True):
    groupName = groupName if groupName is not None else obj.name
    
    
    if type(obj) == list:
        children = []
        for o in obj:
            children += GetParamTreeDict(o) 
    else:
        paramsList = obj.GetParams(asList = True)
        children = []
        for param, value in paramsList:
            if param.guiHide:
                continue
            vType = GetParamType(value) if param.guiType == None else param.guiType
            childrenDict = {"name": param.friendlyName,
                             "type": vType,
                             "value": value,
                             "values": param.guiValues,
                             "obj": obj,
                             "paramName": param.name,
                             "step": param.guiStep,
                             "suffix": param.guiSuffix,
                             "siPrefix": param.guiSiPrefix,
                             "decimals": param.guiDecimals,
                             "limits": param.guiLimits}
            
            for k in childrenDict.keys()[:]:
                if childrenDict[k] is None:
                    del childrenDict[k]
                    
            # print param.name, childrenDict
            
            children.append(childrenDict)
    
    if createGroup:    
        res = [ {'name': groupName, 'type': 'group', 'children': children}]
    else:
        res = children
    return res

        
def GetIndex(arr, value, defaultRes = 0):
    if value in arr:
        return arr.index(value)
    else:
        print "WARNING: Value", value, "not in list"
        return defaultRes


if __name__ == "__main__":
    pass
