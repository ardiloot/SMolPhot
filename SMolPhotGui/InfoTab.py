from PyQt4 import QtGui, QtCore


#===============================================================================
# InfoTab
#===============================================================================

class InfoTab(QtGui.QScrollArea):
    
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
        
        
    def Shown(self):
        self._Update()
        
    def _Update(self):
        
        def AddRow(name, value):
            layout.addWidget(QtGui.QLabel(name), self._UpdateInfoTabRowNr, 0, 1, 1)
            layout.addWidget(QtGui.QLabel(value), self._UpdateInfoTabRowNr, 1, 1, 1)
            self._UpdateInfoTabRowNr += 1
            
        layout = self._layout
        fseries = self.main._fseries
        axialFseries = self.main._axialFseries
        
        # Clear layout
        for i in reversed(range(layout.count())): 
            w = layout.itemAt(i).widget()
            if w is not None:
                w.deleteLater()
            
        self._UpdateInfoTabRowNr = 0    
        if fseries is not None:    
            layout.setColumnStretch(1, 1)
            
            AddRow("Metadata file:", fseries.metadatafile)
            AddRow("Pixel size:", "%.2f nm" % (1e9 * fseries.pixelSize))
            AddRow("", "")
            AddRow("Sequence:", "")
            
            for param, value in fseries.GetParams(asList=True):
                AddRow(param.name, str(value))
            
        if axialFseries is not None:
            AddRow("", "")
            AddRow("Axial calibration:", "")
            
            for param, value in axialFseries.GetParams(asList=True):
                AddRow(param.name, str(value))
            
            spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Minimum, \
                                           QtGui.QSizePolicy.Expanding)
            layout.addItem(spacerItem, self._UpdateInfoTabRowNr, 0, 1, 1)
            

if __name__ == '__main__':
    pass
