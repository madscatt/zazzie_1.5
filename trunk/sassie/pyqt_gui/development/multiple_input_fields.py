__author__ = 'curtisj'
__appname__ = "mif"

import sys
from PyQt4 import QtCore
from PyQt4 import QtGui
import sassie.sasmol.sasmol as sasmol


class myDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(myDialog, self).__init__(parent)

        self.myLabel = QtGui.QLabel("hello")
        self.button = {}
        self.edit = {}
        self.number_of_pairs = 3

        layout = QtGui.QVBoxLayout()
        layout.addWidget(self.myLabel)

        self.dialogOpen()

        myTitle = QtCore.QString("My Title")
        self.setWindowTitle(myTitle)

        for i in xrange(self.number_of_pairs):
            self.button[i] = QtGui.QPushButton()
            st = str(i)
            self.edit[i] = QtGui.QLineEdit(st)
            layout.addWidget(self.edit[i])
            layout.addWidget(self.button[i])

        self.setLayout(layout)


    def dialogOpen(self):
        dialog = Dialog()
        if dialog.exec_():
            self.number_of_pairs = dialog.spinBox.value()
        else:
            QMessageBox.warning(self, __appname__, "Dialog canceled.")

class Dialog(QtGui.QDialog):

    def __init__(self, parent=None):
        super(Dialog, self).__init__(parent)

        self.setWindowTitle("Dialog.")

        self.checkBox = QtGui.QCheckBox("Check me out!")
        self.spinBox = QtGui.QSpinBox()
        buttonOk = QtGui.QPushButton("OK")
        buttonCancel = QtGui.QPushButton("Cancel")

        layout = QtGui.QGridLayout()
        layout.addWidget(self.spinBox, 0, 0)
        layout.addWidget(self.checkBox, 0, 1)
        layout.addWidget(buttonCancel)
        layout.addWidget(buttonOk)
        self.setLayout(layout)

        self.connect(buttonOk, QtCore.SIGNAL("clicked()"), self, QtCore.SLOT("accept()"))
        self.connect(buttonCancel, QtCore.SIGNAL("clicked()"), self, QtCore.SLOT("reject()"))

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    main = myDialog()
    main.show()

    sys.exit(app.exec_())
