from msilib.schema import InstallUISequence
import numpy as np

class Laser():
    def __init__(self, wavelength, k_vector):
        self.wavelength = wavelength
        self.k_vector = k_vector
        
class Instruments():
    def __init__(self):
        pass

class Mirror(Instruments):
    def __init__(self, orientation_vector, reflectivity):
        self.orientation = orientation_vector
        self.reflectivity = reflectivity
    def reflect(self, rayon):
        k_init = rayon.k_vector
        wvl = rayon.wavelength
        mx, my = self.orientation
        kx_init, ky_init = k_init
        if mx != 0:
            theta = np.arctan(my/mx) + np.pi/2
        else :
            theta = 0
        print(theta)
        kx_theta = kx_init*np.cos(2*theta) - ky_init*np.sin(2*theta)
        ky_theta = - ky_init*np.cos(2*theta) - kx_init*np.sin(2*theta)
        k_out = (kx_theta, ky_theta)
        return Laser(wvl, k_out)
    def display(self, position):
        x_position, y_position = position
        x_array = np.linspace(x_position - 1, x_position + 1)
        mx, my = self.orientation
        y_array = (my/mx)*(x_array-x_position) + y_position
        print(y_array)
        return (x_array, y_array)
        # max_length = 2
        # if len(y_array) > 2:
        #     y_array_corrected = y_array[len(y_array)//2]

class TableOptique(dict):
    def __init__(self, size):
        self.size = size

    def add(self, optics, position):
        self[optics] = position

    def draw(self):
        import matplotlib.pyplot as plt
        for optics in self.keys():
            x_array, y_array = optics.display(self[optics])
            plt.plot(x_array, y_array)
        xlim, ylim = self.size
        plt.xlim(-xlim/2, xlim/2)
        plt.ylim(-ylim/2, ylim/2)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()


    




k = (2,3)
rayon = Laser(620, k)
mirror = Mirror((1,1), 1)
table = TableOptique((20,20))
table.add(mirror, (4,4))
table.draw()
# print(mirror.reflect(rayon).k_vector)

