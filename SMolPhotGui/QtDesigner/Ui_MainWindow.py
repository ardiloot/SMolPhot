# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(799, 683)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.centralLayout = QtGui.QGridLayout(self.centralwidget)
        self.centralLayout.setMargin(0)
        self.centralLayout.setObjectName(_fromUtf8("centralLayout"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 799, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuConfig = QtGui.QMenu(self.menubar)
        self.menuConfig.setObjectName(_fromUtf8("menuConfig"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.parametersDockWidget = QtGui.QDockWidget(MainWindow)
        self.parametersDockWidget.setObjectName(_fromUtf8("parametersDockWidget"))
        self.parametersDockWidgetContents = QtGui.QWidget()
        self.parametersDockWidgetContents.setObjectName(_fromUtf8("parametersDockWidgetContents"))
        self.verticalLayout = QtGui.QVBoxLayout(self.parametersDockWidgetContents)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.parametersTree = ParameterTree(self.parametersDockWidgetContents)
        self.parametersTree.setObjectName(_fromUtf8("parametersTree"))
        self.verticalLayout.addWidget(self.parametersTree)
        self.parametersDockWidget.setWidget(self.parametersDockWidgetContents)
        MainWindow.addDockWidget(QtCore.Qt.DockWidgetArea(2), self.parametersDockWidget)
        self.toolBar = QtGui.QToolBar(MainWindow)
        self.toolBar.setIconSize(QtCore.QSize(48, 48))
        self.toolBar.setObjectName(_fromUtf8("toolBar"))
        MainWindow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionOpen = QtGui.QAction(MainWindow)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(_fromUtf8(":/48x48/48x48/Open_48x48.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionOpen.setIcon(icon)
        self.actionOpen.setObjectName(_fromUtf8("actionOpen"))
        self.actionSave = QtGui.QAction(MainWindow)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(_fromUtf8(":/48x48/48x48/Save_48x48.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSave.setIcon(icon1)
        self.actionSave.setObjectName(_fromUtf8("actionSave"))
        self.actionExit = QtGui.QAction(MainWindow)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap(_fromUtf8(":/48x48/48x48/Log Out_48x48.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionExit.setIcon(icon2)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionLoadConfFromFile = QtGui.QAction(MainWindow)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(_fromUtf8(":/48x48/48x48/Stock Index Down_48x48.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionLoadConfFromFile.setIcon(icon3)
        self.actionLoadConfFromFile.setObjectName(_fromUtf8("actionLoadConfFromFile"))
        self.actionSaveConfToFile = QtGui.QAction(MainWindow)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(_fromUtf8(":/48x48/48x48/Stock Index Up_48x48.png")), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSaveConfToFile.setIcon(icon4)
        self.actionSaveConfToFile.setObjectName(_fromUtf8("actionSaveConfToFile"))
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionSave)
        self.menuFile.addAction(self.actionExit)
        self.menuConfig.addSeparator()
        self.menuConfig.addAction(self.actionLoadConfFromFile)
        self.menuConfig.addAction(self.actionSaveConfToFile)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuConfig.menuAction())
        self.toolBar.addAction(self.actionOpen)
        self.toolBar.addAction(self.actionSave)
        self.toolBar.addAction(self.actionLoadConfFromFile)
        self.toolBar.addAction(self.actionSaveConfToFile)
        self.toolBar.addAction(self.actionExit)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "SMolPhot", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuConfig.setTitle(_translate("MainWindow", "Config", None))
        self.parametersDockWidget.setWindowTitle(_translate("MainWindow", "Parameters", None))
        self.toolBar.setWindowTitle(_translate("MainWindow", "toolBar", None))
        self.actionOpen.setText(_translate("MainWindow", "Open", None))
        self.actionSave.setText(_translate("MainWindow", "Save state", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))
        self.actionLoadConfFromFile.setText(_translate("MainWindow", "Load config", None))
        self.actionSaveConfToFile.setText(_translate("MainWindow", "Save config", None))

from pyqtgraph.parametertree import ParameterTree
import resource_rc
