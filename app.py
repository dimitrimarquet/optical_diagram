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

def is_number(n):
    """check if n is a number

    Args:
        n (unknown type)

    Returns:
        bool: True if n is a number, false elsewhise
    """
    try:
        float(n)
        return True
    except ValueError:
        return False

def is_number_pos(n):
    """check if n is a positive number

    Args:
        n (unknown type)

    Returns:
        bool: True if n is a number, false elsewhise
    """
    try:
        float(n)
        if float(n) > 0:
            return True
        return False
    except ValueError:
        return False

    
# from layout_colorwidget import Color
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, width=200, height=200, dpi=100):
        """
        Args:
            width (int, optional): width of the figure. Defaults to 200.
            height (int, optional): height of the figure. Defaults to 200.
            dpi (int, optional): dpi of the figure. Defaults to 100.
        """
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.axes.set_aspect('equal', adjustable='box')
        super().__init__(self.fig)


class MainWindow(QMainWindow):
    def __init__(self):
        
        super().__init__()

        #window parameters
        self.setWindowTitle("Optical diagram")
        self.setFixedWidth(1500)
        self.setFixedHeight(1000)
        
        #layouts of the window
        pagelayout = QVBoxLayout()
        button_layout = QHBoxLayout()
        self.stacklayout = QStackedLayout()
        first_btn_layout = QVBoxLayout()
        zero_btn_layout = QVBoxLayout()
        second_btn_layout = QVBoxLayout()
        third_btn_layout = QVBoxLayout()

        pagelayout.addLayout(button_layout)
        pagelayout.addLayout(self.stacklayout)
        pagelayout.setStretch(0, 3)
        pagelayout.setStretch(1, 7)

        #initialization of the optical problem
        self.table = Optics.TableOptique((20,20))
        self.laser = Optics.Laser(600, (1,1), np.array([1,0,0,0]))
        self.pos_init = (0, 0)

        #set table size
        self.sc = MplCanvas(width=5, height=4, dpi=100)
        self.stacklayout.addWidget(self.sc)

        self.dialog_tab_size = QLineEdit()
        lblName_tab_size = QLabel("Table size (width, height)")
        qbtn_tab_size = QPushButton()
        qbtn_tab_size.setText("Validate table size")
        qbtn_tab_size.clicked.connect(self.set_tab_size)
        lblName_tab_size.setAlignment(Qt.AlignLeft | Qt.AlignBottom)

        zero_btn_layout.addWidget(lblName_tab_size)

        zero_btn_layout.addWidget(self.dialog_tab_size)
        zero_btn_layout.addWidget(qbtn_tab_size)

        zero_btn_layout.setSpacing(0)
        zero_btn_layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        button_layout.addLayout(zero_btn_layout)

        #set laser parameters
        
        #set wavelength
        wvl_btn_layout = QVBoxLayout()
        lblName = QLabel("Wavelength (nm):")
        
        self.dialog = QDoubleSpinBox()
        self.dialog.setMinimum(400)
        self.dialog.setMaximum(700)
        self.dialog.valueChanged.connect(self.set_wvl)

        # Add widgets and correct for alignment
        wvl_btn_layout.addWidget(lblName)
        wvl_btn_layout.addWidget(self.dialog)
        wvl_btn_layout.setContentsMargins(0, 0, 0, 0)
        wvl_btn_layout.setSpacing(0)
        lblName.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        wvl_btn_layout.setAlignment(Qt.AlignTop)
        first_btn_layout.addLayout(wvl_btn_layout)


        # set polarization

        pol_button_layout = QVBoxLayout()
        self.dialog_pol = QLineEdit()
        self.dialog_pol.textEdited.connect(self.enable_validate)
        lblName_pol = QLabel("Polarization (Stokes vector) : ")

        self.qbtn_pol = QPushButton()
        self.qbtn_pol.setText("Validate polarization")
        self.qbtn_pol.clicked.connect(self.set_polarization)

        pol_button_layout.addWidget(lblName_pol)
        pol_button_layout.addWidget(self.dialog_pol)
        pol_button_layout.addWidget(self.qbtn_pol)

        pol_button_layout.setContentsMargins(0, 0, 0, 0)
        pol_button_layout.setSpacing(0)
        lblName_pol.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        pol_button_layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)

        first_btn_layout.addLayout(pol_button_layout)

        #set k vector

        k_button_layout = QVBoxLayout()
        self.dialog_k = QLineEdit()
        self.dialog_k.textEdited.connect(self.enable_validate)

        lblName_k = QLabel("k vector (kx, ky) : ")

        self.qbtn_k = QPushButton()
        self.qbtn_k.setText("Validate k vector")
        self.qbtn_k.clicked.connect(self.set_k_vector)

        k_button_layout.addWidget(lblName_k)
        k_button_layout.addWidget(self.dialog_k)
        k_button_layout.addWidget(self.qbtn_k)
        
        k_button_layout.setContentsMargins(0, 0, 0, 0)
        k_button_layout.setSpacing(0)
        lblName_k.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        first_btn_layout.addLayout(k_button_layout)


        #set source position

        pos_layout = QVBoxLayout()
        self.dialog_pos = QLineEdit()
        self.dialog_pos.textEdited.connect(self.enable_validate)

        lblName_pos = QLabel("Initial position (x, y) : ")

        self.qbtn_pos = QPushButton()
        self.qbtn_pos.setText("Validate position")
        self.qbtn_pos.clicked.connect(self.set_position_init)
        
        pos_layout.addWidget(lblName_pos)

        pos_layout.addWidget(self.dialog_pos)
        pos_layout.addWidget(self.qbtn_pos)
        
        pos_layout.setContentsMargins(0, 0, 0, 0)
        pos_layout.setSpacing(0)
        lblName_pos.setAlignment(Qt.AlignLeft | Qt.AlignBottom)

        first_btn_layout.addLayout(pos_layout)

        #buttons initially disabled to avoid bugs
        self.pol = False
        self.k = False
        self.pos = False

        self.qbtn_k.setEnabled(False)
        self.qbtn_pol.setEnabled(False)
        self.qbtn_pos.setEnabled(False)
        
        first_btn_layout.setSpacing(5)


        #get and set optical elements

        #define a list of optical elements the user can modify
        self.mirror_index = 0
        self.BS_index = 0
        self.polarizer_index = 0

        self.optics_pos = False
        self.optics_orient = False

        self.optics_list = QComboBox()
        self.optics_list.addItems(['Mirror', 'Beam splitter', 'Polarizer'])
        self.optics_list.currentTextChanged.connect(self.enable_add)
        
        #position
        optics_pos = QVBoxLayout()
        self.dialog_pos_optics = QLineEdit()
        self.dialog_pos_optics.textEdited.connect(self.enable_add)
        lblName_pos_optics = QLabel("Position (x, y)")
        lblName_pos_optics.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        optics_pos.setContentsMargins(0, 0, 0, 0)
        optics_pos.setSpacing(0)

        optics_pos.addWidget(lblName_pos_optics)
        optics_pos.addWidget(self.dialog_pos_optics)


        #orientation
        orient_layout = QVBoxLayout()

        self.dialog_orient_optics = QLineEdit()
        self.dialog_orient_optics.textChanged.connect(self.enable_add)
        lblName_orient_optics = QLabel("Orientation of the normal vector (mx, my)")
        # lblName_orient_optics.setText("Orientation")
        lblName_orient_optics.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        orient_layout.setContentsMargins(0, 0, 0, 0)
        orient_layout.setSpacing(0)

        orient_layout.addWidget(lblName_orient_optics)
        orient_layout.addWidget(self.dialog_orient_optics)

        #BS ratio
        BS_ratio_layout = QVBoxLayout()

        self.BS_ratio_dialog = QLineEdit()
        self.BS_ratio_dialog.textChanged.connect(self.enable_add)
        lblName_BS_ratio = QLabel("Beam splitter transmission intensity rate (between 0 and 1)")
        lblName_BS_ratio.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        BS_ratio_layout.setContentsMargins(0, 0, 0, 0)
        BS_ratio_layout.setSpacing(0)

        BS_ratio_layout.addWidget(lblName_BS_ratio)
        BS_ratio_layout.addWidget(self.BS_ratio_dialog)

        self.BS_ratio_dialog.setEnabled(False)
        self.optics_list.currentTextChanged.connect(self.enable_BS)

        #polarizer angle
        pola_angle_layout = QVBoxLayout()

        self.pola_angle_dialog = QLineEdit()
        self.pola_angle_dialog.textChanged.connect(self.enable_add)
        lblName_pola_angle = QLabel("Polarizer angle with respect to horizontal (°)")
        # lblName_orient_optics.setText("Orientation")
        lblName_pola_angle.setAlignment(Qt.AlignLeft | Qt.AlignBottom)
        pola_angle_layout.setContentsMargins(0, 0, 0, 0)
        pola_angle_layout.setSpacing(0)

        pola_angle_layout.addWidget(lblName_pola_angle)
        pola_angle_layout.addWidget(self.pola_angle_dialog)

        self.pola_angle_dialog.setEnabled(False)
        self.optics_list.currentTextChanged.connect(self.enable_pola_angle)

        #buttons to modify the list

        self.btn_add_optics = QPushButton("Add")

        self.btn_add_optics.setEnabled(False)
        self.btn_add_optics.pressed.connect(self.add_optics)

        self.btn_rem_optics = QPushButton("Remove")
        self.btn_rem_optics.pressed.connect(self.remove_optics)

        self.btn_clear_optics = QPushButton("Clear")
        self.btn_clear_optics.pressed.connect(self.clear_optics)
        self.btn_rem_optics.setEnabled(False)
        self.btn_clear_optics.setEnabled(False)

        #add of the layouts in the appropriate order
        second_btn_layout.addWidget(self.optics_list)
        second_btn_layout.addLayout(optics_pos)
        second_btn_layout.addLayout(orient_layout)
        second_btn_layout.addLayout(BS_ratio_layout)
        second_btn_layout.addLayout(pola_angle_layout)

 
        second_btn_layout.addWidget(self.btn_add_optics)

        second_btn_layout.addWidget(self.btn_rem_optics)


        second_btn_layout.addWidget(self.btn_clear_optics)

        button_layout.addLayout(second_btn_layout)

        self.list_widget = QListWidget()

        button_layout.addWidget(self.list_widget)

        button_layout.addLayout(first_btn_layout)


        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)

        #buttons to draw everything and recover the results

        self.btn_draw_optics = QPushButton("Draw Table")
        self.btn_draw_optics.pressed.connect(self.draw_table)
        third_btn_layout.addWidget(self.btn_draw_optics)
        self.btn_draw_optics.setEnabled(False)

        self.btn_draw_laser = QPushButton("Draw Laser")
        self.btn_draw_laser.pressed.connect(self.draw_laser)
        third_btn_layout.addWidget(self.btn_draw_laser)
        self.btn_draw_laser.setEnabled(False)

        btn_clear_table = QPushButton("Clear Table")
        btn_clear_table.pressed.connect(self.clear_table)
        third_btn_layout.addWidget(btn_clear_table)


        btn_dwnl_results = QPushButton("Download")
        btn_dwnl_results.pressed.connect(self.download)
        third_btn_layout.addWidget(btn_dwnl_results)
        button_layout.addLayout(third_btn_layout)



    def set_tab_size(self):
        """set optical table dimensions from user input
        """
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
        elif is_number(text.split()[0]) and  is_number(text.split()[1]):
            x, y = (float(text.split()[0]), float(text.split()[1]))
            self.x_size = [-x/2, x/2]
            self.y_size = [-y/2, y/2]
            self.table.size = (x,y)

            self.sc.axes.set_xlim(self.x_size)
            self.sc.axes.set_ylim(self.y_size)
            self.sc.axes.plot([1], [1])
            self.sc.axes.set_aspect('equal', adjustable='box')

            self.sc.draw()
        self.btn_draw_optics.setEnabled(True)


    def set_wvl(self, i):
        """set wavelength

        Args:
            i (float): wavelength value between 400 and 700 nm (in nm)
        """
        # text = self.dialog.text()
        self.laser.wavelength = i

    def set_polarization(self):
        """set laser polarization
        """
        text = self.dialog_pol.text()
        pol_array = np.array([float(text.split()[0]), float(text.split()[1]), float(text.split()[2]), float(text.split()[3])])
        self.laser.stokes = pol_array
        self.pol = True
        if self.k and self.pos and self.pol :
            self.btn_draw_laser.setEnabled(True)



    def set_k_vector(self):
        """set initial k vector of the laser from user input
        """
        text = self.dialog_k.text()
        tuple = (float(text.split()[0]), float(text.split()[1]))
        if float(text.split()[0]) == 0 and float(text.split()[1]) != 0:
            tuple = (float(text.split()[0])+1e-4, float(text.split()[1]))

        elif float(text.split()[1]) == 0 and float(text.split()[0]) == 0:
            tuple = (float(text.split()[0])+1e-4, float(text.split()[1])+1e-4)
        
        elif float(text.split()[0]) != 0 and float(text.split()[1]) == 0:
            tuple = (float(text.split()[0]), float(text.split()[1])+1e-4)        
        
        self.laser.k_vector = tuple
        self.k = True
        if self.k and self.pos and self.pol :
            self.btn_draw_laser.setEnabled(True)

    def set_position_init(self):
        """set initial position of the laser source from user input
        """
        text = self.dialog_pos.text()
        tuple = (float(text.split()[0]), float(text.split()[1]))
        self.pos_init = tuple
        self.pos = True
        if self.k and self.pos and self.pol :
            self.btn_draw_laser.setEnabled(True)

    def add_optics(self):
        """add optical element to the optical table in the list ['Mirror', 'Beam splitter', 'Polarizer']
        """
        text_pos = self.dialog_pos_optics.text()
        pos_optics = (float(text_pos.split()[0]), float(text_pos.split()[1]))

        text_orient = self.dialog_orient_optics.text()
        orient_optics = (float(text_orient.split()[0]), float(text_orient.split()[1]))

        if float(text_orient.split()[1]) == 0 :
            orient_optics = (float(text_orient.split()[0]), float(text_orient.split()[1])+1e-4)

        if self.optics_list.currentText() == 'Mirror':
            self.table.add(Optics.Mirror(orient_optics, 1, name = 'Mirror ' + str(self.mirror_index)), pos_optics)
            self.list_widget.addItem(self.optics_list.currentText() + str(self.mirror_index))
            self.mirror_index += 1

        elif self.optics_list.currentText() == 'Beam splitter':
            text_ratio = self.BS_ratio_dialog.text()
            ratio = float(text_ratio)

            self.table.add(Optics.BeamSplitter(orient_optics, ratio, name = 'BS' + str(self.BS_index)), pos_optics)
            self.list_widget.addItem(self.optics_list.currentText() + str(self.BS_index))
            self.BS_index += 1

        elif self.optics_list.currentText() == 'Polarizer':
            text_angle = self.pola_angle_dialog.text()
            angle = float(text_angle)
            self.table.add(Optics.Polarizer(orient_optics, angle, name = 'Pol' + str(self.polarizer_index)), pos_optics)
            self.list_widget.addItem(self.optics_list.currentText() + str(self.polarizer_index))
            self.polarizer_index += 1
        
        self.btn_rem_optics.setEnabled(True)
        self.btn_clear_optics.setEnabled(True)

    def enable_add(self):
        """check if every parameter has been entered by the user before adding an optical element and enable or disable the associated button

        Returns:
            None
        """
        txt_orient = self.dialog_orient_optics.text()
        txt_pos = self.dialog_pos_optics.text()

        length_test_orient = len(txt_orient.split()) == 2
        length_test_pos = len(txt_pos.split()) == 2

        if length_test_orient and length_test_pos:
            type_test_orient = is_number(txt_orient.split()[0]) and is_number(txt_orient.split()[1])
            type_test_pos = is_number(txt_pos.split()[0]) and is_number(txt_pos.split()[1])

            if type_test_orient and type_test_pos:

                if self.optics_list.currentText() == 'Mirror':
                    self.btn_add_optics.setEnabled(True)
                    return None

                if self.optics_list.currentText() == 'Beam splitter':
                    text_ratio = self.BS_ratio_dialog.text()
                    if is_number(text_ratio) and 0 <= float(text_ratio) <= 1:
                        self.btn_add_optics.setEnabled(True)
                        return None

                elif self.optics_list.currentText() == 'Polarizer':
                    text_angle = self.pola_angle_dialog.text()
                    if is_number(text_angle) and 0 <= float(text_angle) <= 360 :
                        self.btn_add_optics.setEnabled(True)
                        return None

        self.btn_add_optics.setEnabled(False)
    
    def enable_validate(self):
        """check if every parameter has been entered by the user before validating a laser parameter and enable or disable the associated button
        """
        text_pol = self.dialog_pol.text()
        text_pos_init = self.dialog_pos.text()
        text_k = self.dialog_k.text()
        if len(text_pol.split()) == 4 and is_number_pos(text_pol.split()[0]) and is_number(text_pol.split()[1]) and is_number(text_pol.split()[2]) and is_number(text_pol.split()[3]):
            self.qbtn_pol.setEnabled(True)
        else :
            self.qbtn_pol.setEnabled(False)
        if len(text_pos_init.split()) == 2 and is_number(text_pos_init.split()[0]) and is_number(text_pos_init.split()[1]) :
            self.qbtn_pos.setEnabled(True)
        else :
            self.qbtn_pos.setEnabled(False)

        if len(text_k.split()) == 2 and is_number(text_k.split()[0]) and is_number(text_k.split()[1]) :
            self.qbtn_k.setEnabled(True)
        else :
            self.qbtn_k.setEnabled(False)

    def remove_optics(self):
        """removes the last optical element, disable the button if the new length is 0 and keeps it enabled elsewise
        """
        self.table.popitem()
        row = len(self.list_widget)
        self.list_widget.takeItem(row-1)
        if len(self.table) == 0:
            self.btn_rem_optics.setEnabled(False)
            self.btn_clear_optics.setEnabled(False)

        else:
            self.btn_rem_optics.setEnabled(True)
            self.btn_clear_optics.setEnabled(True)


    def enable_BS(self):
        """check if the optical element the user wants to add is a beam splitter and activates the appropriate parameter inbox
        """
        if self.optics_list.currentText() == 'Beam splitter':
            self.BS_ratio_dialog.setEnabled(True)
        else :
            self.BS_ratio_dialog.setEnabled(False)
        
    def enable_pola_angle(self):
        """check if the optical element the user wants to add is a polarizer and activates the appropriate parameter inbox
        """
        if self.optics_list.currentText() == 'Polarizer':
            self.pola_angle_dialog.setEnabled(True)
        else :
            self.pola_angle_dialog.setEnabled(False)


    def clear_optics(self):
        """removes everything from the optical element list
        """
        self.list_widget.clear()
        self.table.clear()
        self.btn_rem_optics.setEnabled(False)
        self.btn_clear_optics.setEnabled(False)



    def draw_table(self):
        """draws the optical table
        """
        x_lim = self.sc.axes.get_xlim()
        y_lim = self.sc.axes.get_ylim()

        self.sc.axes.clear()
        self.table.draw(self.sc.axes)

        self.sc.axes.set_xlim(x_lim)
        self.sc.axes.set_ylim(y_lim)
        self.sc.axes.set_aspect('equal', adjustable='box')

        self.sc.draw()

    def draw_laser(self):
        """draws the laser path
        """
        x_lim = self.sc.axes.get_xlim()
        y_lim = self.sc.axes.get_ylim()

        # self.sc.axes.clear()
        self.table.draw_laser(self.laser,self.pos_init, self.sc.axes)

        self.sc.axes.set_xlim(x_lim)
        self.sc.axes.set_ylim(y_lim)
        self.sc.axes.set_aspect('equal', adjustable='box')

        self.sc.draw()
    
    def clear_table(self):
        """removes everything drawn on the optical table
        """
        x_lim = self.sc.axes.get_xlim()
        y_lim = self.sc.axes.get_ylim()
        self.sc.axes.clear()

        self.sc.axes.set_xlim(x_lim)
        self.sc.axes.set_ylim(y_lim)
        self.sc.axes.set_aspect('equal', adjustable='box')

        self.sc.draw()

    def download(self):
        """downloads the results : png of the drawing and report.txt file
        """
        self.table.report(self.laser, self.pos_init, filename = "rapport_optique.txt")
        self.sc.fig.savefig("optical_path.jpg")


