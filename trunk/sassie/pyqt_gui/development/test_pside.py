__author__ = 'curtisj'
import sys

import PySide
#from PySide.QtGui import *

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import numpy

# Create the application object

app = QApplication(sys.argv)

# Create a simple dialog box
msgBox = QMessageBox()
text = 'hey this works! ' + __author__ + ' '+str(numpy.sin((3.2222)))
msgBox.setText(text + PySide.__version__)


msgBox.exec_()


