import sys
import Optics
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from PyQt5.QtWidgets import (
    QApplication,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QStackedLayout,
    QVBoxLayout,
    QWidget,
    QListWidget,
)

# from layout_colorwidget import Color
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=200, height=200, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        self.axes.set_xlim(-3,4)
        self.axes.set_ylim(-3,4)

        super().__init__(fig)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Optical diagram")
        self.setFixedWidth(1000)
        self.setFixedHeight(1000)

        pagelayout = QVBoxLayout()
        button_layout = QHBoxLayout()
        self.stacklayout = QStackedLayout()

        first_btn_layout = QVBoxLayout()

        pagelayout.addLayout(button_layout)
        pagelayout.addLayout(self.stacklayout)

        self.sc = MplCanvas(self, width=5, height=4, dpi=100)

        self.xinit = [0,1]
        self.yinit = [0,1]
        self.x = self.xinit
        self.y = self.yinit
        # self.sc.axes.plot(self.x, self.y)

        # btn = QPushButton("red")
        # btn.pressed.connect(self.activate_tab_1)
        # button_layout.addWidget(btn)
        self.stacklayout.addWidget(self.sc)

        btn = QPushButton("red_bis")
        btn.pressed.connect(self.activate_tab_1)
        first_btn_layout.addWidget(btn)

        btn = QPushButton("red_bat")
        btn.pressed.connect(self.activate_tab_1)
        first_btn_layout.addWidget(btn)

        button_layout.addLayout(first_btn_layout)

        self.btn2 = QPushButton("green")
        self.btn2.pressed.connect(self.activate_tab_2)
        button_layout.addWidget(self.btn2)
        # self.stacklayout.addWidget(Color("green"))

        # if btn.pressed.emit():
        #     print('tret')
        #     btn2.setEnabled(False)
        self.list_widget = QListWidget()
        # self.list_widget.resize(100,100)
        button_layout.addWidget(self.list_widget)

        btn = QPushButton("yellow")
        btn.pressed.connect(self.activate_tab_3)
        button_layout.addWidget(btn)
        # self.stacklayout.addWidget(Color("yellow"))
        

        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)


    def activate_tab_1(self):
        print('test1')
        print(self.setWindowTitle('test'))
        self.sc.axes.clear()
        table = Optics.TableOptique((20,20))
        self.sc.axes.set_xlim(-table.size[0]/2, table.size[0]/2)
        self.sc.axes.set_ylim(-table.size[1]/2, table.size[1]/2)

        table.add(Optics.Mirror((1,1), 1), (0,0))
        print('r')
        table.draw(self.sc.axes)
        self.sc.draw()
        self.btn2.setEnabled(False)
        self.list_widget.addItem(Optics.Mirror((1,1), 1).name)
        return True

    def activate_tab_2(self):
        self.stacklayout.setCurrentIndex(1)

    def activate_tab_3(self):
        self.stacklayout.setCurrentIndex(2)


app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()