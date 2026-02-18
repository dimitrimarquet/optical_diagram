from PyQt5.QtWidgets import QApplication, QWidget, QHBoxLayout, QPushButton, QSizePolicy

app = QApplication([])

window = QWidget()
layout = QHBoxLayout()

# Add widgets
button1 = QPushButton("Button 1")
button2 = QPushButton("Button 2")
button3 = QPushButton("Button 3")

# Allow widgets to expand
button1.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
button2.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
button3.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

layout.addWidget(button1)
layout.addWidget(button2)
layout.addWidget(button3)

# Set stretch factors
layout.setStretch(0, 1)
layout.setStretch(1, 2)
layout.setStretch(2, 46)

window.setLayout(layout)
window.show()
app.exec_()