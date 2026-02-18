import app
import sys

appli = app.QApplication(sys.argv)

window = app.MainWindow()   
window.show()

appli.exec()