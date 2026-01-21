from msilib.schema import InstallUISequence
import numpy as np
import matplotlib.pyplot as plt

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
        kx_theta = kx_init*np.cos(2*theta) - ky_init*np.sin(2*theta)
        ky_theta = - ky_init*np.cos(2*theta) - kx_init*np.sin(2*theta)
        k_out = (kx_theta, ky_theta)
        print(k_out)
        return Laser(wvl, k_out)
   
    # def display(self, position):
    #     x_position, y_position = position
    #     x_array = np.linspace(x_position - 1, x_position + 1)
    #     mx, my = self.orientation
    #     if mx != 0:
    #         y_array = (my/mx)*(x_array-x_position) + y_position
    #         return (x_array, y_array)
    #     else :
    #         y_array = np.linspace(y_position - 1, y_position +1) 
    #         return (x_array, y_array)
        # print(y_array)
        # max_length = 2
        # if len(y_array) > 2:
        #     y_array_corrected = y_array[len(y_array)//2]
    def display(self, position, length = 2):
        x_position, y_position = position
        mx, my = self.orientation
        if mx != 0:
            theta = np.arctan(my/mx) + np.pi/2
        else :
            theta = 0
        x_max = x_position + np.cos(theta)*length/2
        x_min = x_position - np.cos(theta)*length/2

        y_max = x_position + np.sin(theta)*length/2
        y_min = x_position - np.sin(theta)*length/2

        x_array = np.linspace(x_min, x_max, 10)
        y_array = np.linspace(y_min, y_max, 10)
        return (x_array, y_array)

class TableOptique(dict):
    def __init__(self, size):
        self.size = size

    def add(self, optics, position):
        self[optics] = position

    def draw(self):
        for optics in self.keys():
            x_array, y_array = optics.display(self[optics])
            plt.scatter(x_array, y_array)
        xlim, ylim = self.size
        plt.xlim(-xlim/2, xlim/2)
        plt.ylim(-ylim/2, ylim/2)
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()
    
    def path_laser(self, laser, position_init, epsilon = 0.005):
        ray = laser
        x_init, y_init = position_init
        x_size, y_size = self.size
        x_array = [x_init]
        y_array = [y_init]
        d_pos = 0.01
        
        while x_array[-1] < x_size/2 and y_array[-1] < y_size/2:
            kx_init, ky_init = ray.k_vector
            x_array.append(x_array[-1]+kx_init*d_pos)
            y_array.append(y_array[-1]+ky_init*d_pos)
            for optics in self.keys():
                position_optics = self[optics]
                x_optics, y_optics = optics.display(position_optics)
                for x in x_optics:
                    for y in y_optics :
                        # print(x_array[-1], x)
                        if np.sqrt((x_array[-1] - x)**2 + (y_array[-1] - y)**2) < epsilon :
                        # if x_array[-1] > x - epsilon and x_array[-1] < x + epsilon :
                        #     if y_array[-1] > y - epsilon and x_array[-1] < y + epsilon:
                            print(x_array[-1], x)

                            ray = optics.reflect(ray)
                        

        return (x_array, y_array)
    
    def draw_laser(self, laser, position_source):
        x_array, y_array = self.path_laser(laser, position_source)
        plt.scatter(x_array, y_array)
        # plt.show()

k = (2,1)
rayon = Laser(620, k)
mirror1 = Mirror((2,3), 1)
# mirror2 = Mirror((1,1), 1)
table = TableOptique((20,20))
table.add(mirror1, (4,4))
# table.add(mirror2, (4,4))

table.draw()
table.draw_laser(rayon, (-1,1))
# print(mirror.reflect(rayon).k_vector)

plt.show()