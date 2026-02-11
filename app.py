import sys
import Optics
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import numpy as np
from PyQt5.QtCore import Qt

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
    QLineEdit,
    QSpinBox,
    QComboBox,
    QDoubleSpinBox,
    QGroupBox
)

# from layout_colorwidget import Color
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=200, height=200, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.axes.set_aspect('equal', adjustable='box')
        # self.axes.set_xlim(-3,4)
        # self.axes.set_ylim(-3,4)

        super().__init__(self.fig)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Optical diagram")
        self.setFixedWidth(1500)
        self.setFixedHeight(1000)
        self.mirror_index = 0
        self.BS_index = 0
        self.polarizer_index = 0

        pagelayout = QVBoxLayout()
        button_layout = QHBoxLayout()
        self.stacklayout = QStackedLayout()
        first_btn_layout = QVBoxLayout()
        zero_btn_layout = QVBoxLayout()
        second_btn_layout = QVBoxLayout()
        third_btn_layout = QVBoxLayout()

        pagelayout.addLayout(button_layout)
        pagelayout.addLayout(self.stacklayout)

        self.table = Optics.TableOptique((20,20))
        self.laser = Optics.Laser(600, (1,1), np.array([1,0,0,0]))
        self.pos_init = (0, 0)

        #set table size

        self.sc = MplCanvas(self, width=5, height=4, dpi=100)
        self.stacklayout.addWidget(self.sc)

        self.dialog_tab_size = QLineEdit()
        lblName_tab_size = QLabel(self.dialog_tab_size)
        lblName_tab_size.setText("Size table : ")
        qbtn_tab_size = QPushButton()
        qbtn_tab_size.setText("Validate")
        qbtn_tab_size.clicked.connect(self.set_tab_size)

        zero_btn_layout.addWidget(self.dialog_tab_size)
        zero_btn_layout.addWidget(qbtn_tab_size)
        button_layout.addLayout(zero_btn_layout)

        #get wvl
        # self.dialog = QLineEdit()
        # lblName = QLabel(self.dialog)
        # lblName.setText("Wavelength : ")
        # qbtn = QPushButton()
        # qbtn.setText("Validate")
        # qbtn.clicked.connect(self.set_wvl)

        wvl_btn_layout = QVBoxLayout()
        
        lblName = QLabel("Wavelength (nm):")
        self.dialog = QDoubleSpinBox()

        self.dialog.setMinimum(400)
        self.dialog.setMaximum(700)

        self.dialog.valueChanged.connect(self.set_wvl)

        # Add widgets
        wvl_btn_layout.addWidget(lblName)
        wvl_btn_layout.addWidget(self.dialog)

        # 🔥 Remove all extra space
        wvl_btn_layout.setContentsMargins(0, 0, 0, 0)
        wvl_btn_layout.setSpacing(0)

        # Align label tightly
        lblName.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        # self.dialog.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        first_btn_layout.addLayout(wvl_btn_layout)

        # first_btn_layout.addWidget(qbtn)

        # get polarization

        pol_button_layout = QVBoxLayout()
        self.dialog_pol = QLineEdit()
        lblName_pol = QLabel("Polarization : ")
        # lblName_pol.setText("Polarization : ")

        qbtn_pol = QPushButton()
        qbtn_pol.setText("Validate polarization")
        qbtn_pol.clicked.connect(self.set_polarization)

        # pol_button_layout.addWidget()
        pol_button_layout.addWidget(lblName_pol)
        pol_button_layout.addWidget(self.dialog_pol)
        pol_button_layout.addWidget(qbtn_pol)

        # 🔥 Remove all extra space
        pol_button_layout.setContentsMargins(0, 0, 0, 0)
        pol_button_layout.setSpacing(0)
        lblName_pol.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        first_btn_layout.addLayout(pol_button_layout)

        #get k vector

        # self.dialog_k = QLineEdit()
        # lblName_k = QLabel(self.dialog_k)
        # lblName_k.setText("k-vector : ")
        # qbtn_k = QPushButton()
        # qbtn_k.setText("Validate")
        # qbtn_k.clicked.connect(self.set_k_vector)

        # # self.dialog_k = QDoubleSpinBox()
        # # lblName_k = QLabel(self.)
        # # lblName_k.setText("Wavelength : ")

        # first_btn_layout.addWidget(self.dialog_k)
        # first_btn_layout.addWidget(qbtn_k)

        k_button_layout = QVBoxLayout()
        self.dialog_k = QLineEdit()
        lblName_k = QLabel("k vector : ")
        # lblName_pol.setText("Polarization : ")

        qbtn_k = QPushButton()
        qbtn_k.setText("Validate polarization")
        qbtn_k.clicked.connect(self.set_k_vector)

        # pol_button_layout.addWidget()
        k_button_layout.addWidget(lblName_k)
        k_button_layout.addWidget(self.dialog_k)
        k_button_layout.addWidget(qbtn_k)
        
        # 🔥 Remove all extra space
        k_button_layout.setContentsMargins(0, 0, 0, 0)
        k_button_layout.setSpacing(0)
        lblName_k.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        first_btn_layout.addLayout(k_button_layout)


        #get initial position

        pos_layout = QVBoxLayout()
        self.dialog_pos = QLineEdit()
        lblName_pos = QLabel("position_init : ")
        # lblName_pos.setText("position_init : ")

        qbtn_pos = QPushButton()
        qbtn_pos.setText("Validate position")
        qbtn_pos.clicked.connect(self.set_position_init)
        
        pos_layout.addWidget(lblName_pos)

        pos_layout.addWidget(self.dialog_pos)
        pos_layout.addWidget(qbtn_pos)
        
        pos_layout.setContentsMargins(0, 0, 0, 0)
        pos_layout.setSpacing(0)
        lblName_pos.setAlignment(Qt.AlignLeft | Qt.AlignBottom)

        first_btn_layout.addLayout(pos_layout)
        button_layout.addLayout(first_btn_layout)

        #get and set optical elements

        self.list_widget = QListWidget()
        # self.list_widget.resize(100,100)
        button_layout.addWidget(self.list_widget)

        self.optics_list = QComboBox()
        self.optics_list.addItems(['Mirror', 'Beam splitter', 'Polarizer'])
        
        self.dialog_pos_optics = QLineEdit()
        lblName_pos_optics = QLabel(self.dialog_pos_optics)
        lblName_pos_optics.setText("Position optics")

        self.dialog_orient_optics = QLineEdit()
        lblName_orient_optics = QLabel(self.dialog_orient_optics)
        lblName_orient_optics.setText("Orientation")

        btn_add_optics = QPushButton("Add")
        btn_add_optics.pressed.connect(self.add_optics)

        # btn_add_optics = QPushButton("Remove")
        # btn_add_optics.pressed.connect(self.add_optics)

        btn_clear_optics = QPushButton("Clear")
        btn_clear_optics.pressed.connect(self.clear_optics)
        
        second_btn_layout.addWidget(self.optics_list)
        second_btn_layout.addWidget(self.dialog_pos_optics)
        second_btn_layout.addWidget(self.dialog_orient_optics)

 
        second_btn_layout.addWidget(btn_add_optics)
        second_btn_layout.addWidget(btn_clear_optics)

        button_layout.addLayout(second_btn_layout)


        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)

        btn_draw_optics = QPushButton("Draw Table")
        btn_draw_optics.pressed.connect(self.draw_table)
        third_btn_layout.addWidget(btn_draw_optics)

        btn_draw_laser = QPushButton("Draw Laser")
        btn_draw_laser.pressed.connect(self.draw_laser)
        third_btn_layout.addWidget(btn_draw_laser)

        btn_dwnl_results = QPushButton("Download")
        btn_dwnl_results.pressed.connect(self.download)
        third_btn_layout.addWidget(btn_dwnl_results)
        button_layout.addLayout(third_btn_layout)

        #drawing the laser path


    def set_tab_size(self):
        text = self.dialog_tab_size.text()
        self.sc.axes.clear()
        if len(text) == 0:
            self.sc.axes.clear()
            self.x_size = [0,1]
            self.y_size = [0,1]
            self.sc.axes.set_xlim(self.x_size)
            self.sc.axes.set_ylim(self.y_size)

            self.sc.axes.plot([1], [1])
            self.sc.axes.set_aspect('equal', adjustable='box')
            self.sc.draw()

            self.table.size = (1,1)
        else :
            x, y = (float(text.split()[0]), float(text.split()[1]))
            self.x_size = [-x/2, x/2]
            self.y_size = [-y/2, y/2]
            self.table.size = (x,y)

            self.sc.axes.set_xlim(self.x_size)
            self.sc.axes.set_ylim(self.y_size)
            self.sc.axes.plot([1], [1])
            self.sc.axes.set_aspect('equal', adjustable='box')

            self.sc.draw()


    def set_wvl(self, i):
        # text = self.dialog.text()
        self.laser.wavelength = i

    def set_polarization(self):
        text = self.dialog_pol.text()
        pol_array = np.array([float(text.split()[0]), float(text.split()[1]), float(text.split()[2]), float(text.split()[3])])
        self.laser.stokes_vector = pol_array


    def set_k_vector(self):
        text = self.dialog_k.text()
        tuple = (float(text.split()[0]), float(text.split()[1]))
        self.laser.k_vector = tuple

    def set_position_init(self):
        text = self.dialog_pos.text()
        tuple = (float(text.split()[0]), float(text.split()[1]))
        self.pos_init = tuple

    # def set_k_vector(self, k_vector):


    def add_optics(self):
        # print('test1')
        # print(self.setWindowTitle('test'))
        # self.sc.axes.clear()
        # table = Optics.TableOptique((20,20))
        # self.sc.axes.set_xlim(-table.size[0]/2, table.size[0]/2)
        # self.sc.axes.set_ylim(-table.size[1]/2, table.size[1]/2)

        # table.add(Optics.Mirror((1,1), 1), (0,0))
        # print('r')
        # table.draw(self.sc.axes)
        # self.sc.draw()
        # self.btn2.setEnabled(False)
        # self.list_widget.addItem(Optics.Mirror((1,1), 1).name)
        # return True
        
        text_pos = self.dialog_pos_optics.text()
        pos_optics = (float(text_pos.split()[0]), float(text_pos.split()[1]))

        text_orient = self.dialog_orient_optics.text()
        orient_optics = (float(text_orient.split()[0]), float(text_orient.split()[1]))
        if self.optics_list.currentText() == 'Mirror':
            self.table.add(Optics.Mirror(orient_optics, 1, name = 'Mirror ' + str(self.mirror_index)), pos_optics)
            self.list_widget.addItem(self.optics_list.currentText() + str(self.mirror_index))
            self.mirror_index += 1

        elif self.optics_list.currentText() == 'Beam splitter':
            self.table.add(Optics.BeamSplitter(orient_optics, 0.5), pos_optics)
            self.list_widget.addItem(self.optics_list.currentText() + str(self.BS_index))
            self.BS_index += 1

        elif self.optics_list.currentText() == 'Polarizer':
            self.table.add(Optics.Polarizer(orient_optics, 1), pos_optics)
            self.list_widget.addItem(self.optics_list.currentText() + str(self.polarizer_index))
            self.polarizer_index += 1


    # def remove_optics(self):


    def clear_optics(self):
        self.list_widget.clear()
        self.table.clear()

    def draw_table(self):
        x_lim = self.sc.axes.get_xlim()
        y_lim = self.sc.axes.get_ylim()

        self.sc.axes.clear()
        self.table.draw(self.sc.axes)

        self.sc.axes.set_xlim(x_lim)
        self.sc.axes.set_ylim(y_lim)
        self.sc.axes.set_aspect('equal', adjustable='box')

        self.sc.draw()

    def draw_laser(self):
        x_lim = self.sc.axes.get_xlim()
        y_lim = self.sc.axes.get_ylim()

        # self.sc.axes.clear()
        self.table.draw_laser(self.laser,self.pos_init, self.sc.axes)

        self.sc.axes.set_xlim(x_lim)
        self.sc.axes.set_ylim(y_lim)
        self.sc.axes.set_aspect('equal', adjustable='box')

        self.sc.draw()
    
    def download(self):
        self.table.report(self.laser, self.pos_init, filename = "rapport_optique.tx")
        self.sc.fig.savefig("optical_path.jpg")


app = QApplication(sys.argv)

window = MainWindow()   
window.show()

app.exec()