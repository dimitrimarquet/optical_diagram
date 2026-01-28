import sys
import sys
from PyQt5.QtWidgets import QWidget, QPushButton, QVBoxLayout, QApplication, QMainWindow, QSlider, QDoubleSpinBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np

def Mrot(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        self.axes.set_xlim(-3,4)
        self.axes.set_ylim(-3,4)

        super().__init__(fig)

class Mirror():
    def __init__(self, x, y):
        self.pos = (x, y)
    def draw(self):
        x, y = self.pos
        x_array, y_array = ([x-1/2, x + 1/2], [y - 1/2, y + 1/2])
        return x_array, y_array

class MainWindow(QMainWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Create the maptlotlib FigureCanvas object,
        # which defines a single set of axes as self.axes.
        self.sc = MplCanvas(self, width=5, height=4, dpi=100)

        self.xinit = [0,1]
        self.yinit = [0,1]
        self.x = self.xinit
        self.y = self.yinit
        self.sc.axes.plot(self.x, self.y)
        # self.button_print_1 = QPushButton('Bouton surprise 1', self)
        # self.button_print_1.setFixedSize(150, 30)
        # self.texte = 'Bouton surprise 1'
        # self.button_print_1.clicked.connect(self.print)
        layout = QVBoxLayout()
        
        widget = QSlider()
        layout.setSpacing(20)
        widget.setMinimum(-200)
        widget.setMaximum(200)
        widget.setSingleStep(1)
        widget.valueChanged.connect(self.value_changed)
        layout.addWidget(widget)

        button = QPushButton('Bouton surprise 1', self)
        self.optics = []
        button.clicked.connect(self.draw_mirror)
        layout.addWidget(button)

        layout.addWidget(self.sc)

        spinBox = QDoubleSpinBox()
        spinBox.setRange(0,2*np.pi)
        spinBox.setSingleStep(0.1)
        spinBox.valueChanged.connect(self.change_direction)
        layout.addWidget(spinBox)

        widget = QWidget()
        widget.setLayout(layout)

        self.setCentralWidget(widget)

    def value_changed(self, i):
        self.sc.axes.clear()
        # self.setCentralWidget(self.sc)
        for o in self.optics :
            x, y = o.draw()
            self.sc.axes.plot(x, y)

        norm = np.sqrt((self.x[1]-self.x[0])**2 + (self.y[1]-self.y[0])**2)
        norm_init = np.sqrt((self.xinit[1]-self.xinit[0])**2 + (self.yinit[1]-self.yinit[0])**2)
        print(norm, norm_init)
        # x = np.array([self.x, self.y])/norm_init
        # y = np.array([self.x, self.y])/norm_init
        rap = self.y[1]/self.x[1]
        self.x = [self.x[0],i*self.x[1]*norm_init/norm*1/200]
        self.y = [self.y[0],i*self.y[1]*norm_init/norm*rap/200]

        self.sc.axes.plot(self.x, self.y)
        self.sc.axes.set_xlim(-3,4)
        self.sc.axes.set_ylim(-3,4)

        self.sc.draw()

        print("checj", i*self.x[1]/200)

        print(i)
    
    def draw_mirror(self):
        self.sc.axes.clear()
        m1 = Mirror(1,1)
        self.optics.append(m1)
        x, y = m1.draw()
        self.sc.axes.plot(x, y)
        self.sc.axes.set_xlim(-3,4)
        self.sc.axes.set_ylim(-3,4)
        self.sc.axes.plot(self.x, self.y)

        self.sc.draw()

    def change_direction(self, theta):
        self.sc.axes.clear()
        # self.setCentralWidget(self.sc)
        for o in self.optics :
            x, y = o.draw()
            self.sc.axes.plot(x, y)

        norm = np.sqrt((self.x[1]-self.x[0])**2 + (self.y[1]-self.y[0])**2)
        norm_init = np.sqrt((self.xinit[1]-self.xinit[0])**2 + (self.yinit[1]-self.yinit[0])**2)
        vect_rot = Mrot(theta) @ np.array([self.xinit[1]-self.xinit[0], self.yinit[1]-self.yinit[0]])*norm/norm_init
        self.x = [self.xinit[0], vect_rot[0] + self.xinit[0]]
        self.y = [self.yinit[0], vect_rot[1] + self.y[0]]
        print(theta, self.x)
        self.sc.axes.plot(self.x, self.y)
        self.sc.axes.set_xlim(-3,4)
        self.sc.axes.set_ylim(-3,4)

        self.sc.draw()

        print("checj")




m1 = Mirror(1,1) 
app = QApplication(sys.argv)
w = MainWindow()
w.show()
app.exec_()

# fig, ax = plt.subplots(121)
# plt.plot([0,1], [0,1])
# plt.show()