from msilib.schema import InstallUISequence
import numpy as np
import matplotlib.pyplot as plt

class Laser():
    """
    Représente un faisceau laser avec une longueur d'onde et un vecteur d'onde.

    Attributs:
        wavelength (float): Longueur d'onde du laser (en nm ou m selon l'unité choisie).
        k_vector (tuple): Vecteur d'onde (k_x, k_y) représentant la direction de propagation.
    """
    def __init__(self, wavelength, k_vector):
        self.wavelength = wavelength
        self.k_vector = k_vector

class Instruments():
    """
    Classe de base pour tous les instruments optiques.
    """
    def __init__(self):
        pass

class Mirror(Instruments):
    """
    Représente un miroir optique avec une orientation et une réflectivité.

    Attributs:
        orientation (tuple): Vecteur d'orientation du miroir (mx, my).
        reflectivity (float): Coefficient de réflectivité du miroir (0 à 1).
    """
    def __init__(self, orientation_vector, reflectivity):
        self.orientation = orientation_vector
        self.reflectivity = reflectivity

    def reflect(self, rayon):
        """
        Réfléchit un faisceau laser incident selon la loi de la réflexion.

        Args:
            rayon (Laser): Faisceau laser incident.

        Returns:
            Laser: Faisceau laser réfléchi avec un nouveau vecteur d'onde.
        """
        k_init = rayon.k_vector
        wvl = rayon.wavelength
        mx, my = self.orientation
        kx_init, ky_init = k_init
        if mx != 0:
            theta = -np.pi/2 + np.arctan(my/mx)
        else :
            theta = 0
        kx_theta = kx_init*np.cos(2*theta) + ky_init*np.sin(2*theta)
        ky_theta = - ky_init*np.cos(2*theta) + kx_init*np.sin(2*theta)
        k_out = (kx_theta, ky_theta)
        return Laser(wvl, k_out)

    def display(self, position, length = 2):
        """
        Calcule les coordonnées pour afficher le miroir sous forme de segment.

        Args:
            position (tuple): Position (x, y) du centre du miroir.
            length (float): Longueur du miroir à afficher.

        Returns:
            tuple: Deux tableaux numpy (x_array, y_array) représentant les extrémités du miroir.
        """
        x_position, y_position = position
        mx, my = self.orientation
        if mx != 0:
            theta = - np.pi/2 + np.arctan(my/mx)
        else :
            theta = 0
        x_max = x_position + np.cos(theta)*length/2
        x_min = x_position - np.cos(theta)*length/2

        y_max = y_position + np.sin(theta)*length/2
        y_min = y_position - np.sin(theta)*length/2

        x_array = np.linspace(x_min, x_max, 10)
        y_array = np.linspace(y_min, y_max, 10)
        return (x_array, y_array)

class TableOptique(dict):
    """
    Représente une table optique contenant divers instruments optiques.

    Attributs:
        size (tuple): Taille (largeur, hauteur) de la table.
    """
    def __init__(self, size):
        self.size = size

    def add(self, optics, position):
        """
        Ajoute un instrument optique à la table.

        Args:
            optics (Instruments): Instrument optique à ajouter.
            position (tuple): Position (x, y) de l'instrument sur la table.
        """
        self[optics] = position

    def draw(self):
        """
        Affiche tous les instruments optiques présents sur la table.
        """
        for optics in self.keys():
            x_array, y_array = optics.display(self[optics])
            plt.plot(x_array, y_array, color = "red")
        xlim, ylim = self.size
        plt.xlim(-xlim/2, xlim/2)
        plt.ylim(-ylim/2, ylim/2)
        plt.gca().set_aspect('equal', adjustable='box')
        # plt.show()

    def path_laser(self, laser, position_init, epsilon = 0.005):
        """
        Calcule le trajet d'un faisceau laser sur la table optique.

        Args:
            laser (Laser): Faisceau laser à propager.
            position_init (tuple): Position initiale (x, y) du faisceau.
            epsilon (float): Précision pour la détection des intersections.

        Returns:
            tuple: Deux listes (x_array, y_array) représentant le trajet du faisceau.
        """
        ray = laser
        x_init, y_init = position_init
        x_size, y_size = self.size
        x_array = [x_init]
        y_array = [y_init]
        d_pos = 15
        used_optics = []
        while np.abs(x_array[-1]) < x_size/2 and np.abs(y_array[-1]) < y_size/2 :
            list_optics_in_cone = []
            kx_init, ky_init = ray.k_vector
            norm_init = np.sqrt(kx_init**2 + ky_init**2)

            for optics in self.keys():
                x_optics, y_optics = self[optics]
                x_array_optics, y_array_optics = optics.display(self[optics])

                i_xmax_optics = np.argmax(x_array_optics)
                i_xmin_optics = np.argmin(x_array_optics)
                i_ymax_optics = np.argmax(y_array_optics)
                i_ymin_optics = np.argmin(y_array_optics)

                kmax_x, kmax_y = (x_array_optics[i_xmax_optics] - x_init, y_array_optics[i_ymax_optics] - y_init)
                kmin_x, kmin_y = (x_array_optics[i_xmin_optics] - x_init, y_array_optics[i_ymin_optics] - y_init)

                norm_max = np.sqrt(kmax_x**2 + kmax_y**2)
                norm_min = np.sqrt(kmin_x**2 + kmin_y**2)
                testy = kmin_y/norm_min >= ky_init/norm_init >= kmax_y/norm_max or kmax_y/norm_max >= ky_init/norm_init >= kmin_y/norm_min
                test_optics_in_cone = kmin_x/norm_min >= kx_init/norm_init >= kmax_x/norm_max or kmin_x/norm_min <= kx_init/norm_init <= kmax_x/norm_max

                if test_optics_in_cone and testy and optics not in used_optics:
                    list_optics_in_cone.append(optics)
            if len(list_optics_in_cone) != 0:
                list_position_optics_in_cone = [self[el] for el in list_optics_in_cone]
                print(list_position_optics_in_cone)
                list_distance_optics_in_cone = [np.sqrt((el[0] - x_array[-1])**2 + (el[1] - y_array[-1])**2) for el in list_position_optics_in_cone]

                optics = list_optics_in_cone[np.argmin(list_distance_optics_in_cone)]
                position_optics = self[optics]
                x_optics, y_optics = position_optics
                x_array_optics, y_array_optics = optics.display(position_optics)

                i_xmax_optics = np.argmax(x_array_optics)
                i_xmin_optics = np.argmin(x_array_optics)
                i_ymax_optics = np.argmax(y_array_optics)
                i_ymin_optics = np.argmin(y_array_optics)

                kmax_x, kmax_y = (x_array_optics[i_xmax_optics] - x_init, y_array_optics[i_ymax_optics] - y_init)
                kmin_x, kmin_y = (x_array_optics[i_xmin_optics] - x_init, y_array_optics[i_ymin_optics] - y_init)

                norm_max = np.sqrt(kmax_x**2 + kmax_y**2)
                norm_min = np.sqrt(kmin_x**2 + kmin_y**2)

                if x_array_optics[i_xmax_optics] - x_array_optics[i_xmin_optics] != 0:
                    coef_miroir = (y_array_optics[i_xmax_optics] - y_array_optics[i_xmin_optics])/(x_array_optics[i_xmax_optics] - x_array_optics[i_xmin_optics])
                    coef_laser = ky_init/kx_init

                    if coef_miroir-coef_laser != 0:
                        x_intersect = (y_init-y_optics + coef_miroir*x_optics - coef_laser*x_init)/(coef_miroir-coef_laser)
                        y_intersect = coef_miroir*(x_intersect-x_optics) + y_optics
                        x_array.append(x_intersect)
                        y_array.append(y_intersect)
                        x_init, y_init = x_intersect, y_intersect
                        ray = optics.reflect(ray)
                        used_optics.append(optics)

                    else :
                        x_array.append(x_optics)
                        y_array.append(y_optics)
                        return (x_array, y_array)
                else :
                    used_optics.append(optics)
                    print("not implemented yet...")

            else:
                x_array.append(x_array[-1] + d_pos*kx_init)
                y_array.append(y_array[-1] + d_pos*ky_init)

        return (x_array, y_array)

    def draw_laser(self, laser, position_source):
        """
        Affiche le trajet d'un faisceau laser sur la table optique.

        Args:
            laser (Laser): Faisceau laser à afficher.
            position_source (tuple): Position initiale (x, y) du faisceau.
        """
        x_array, y_array = self.path_laser(laser, position_source)
        plt.plot(x_array, y_array)
        # plt.show()

# Exemple d'utilisation
k = (-5,3.5)
rayon = Laser(620, k)
mirror1 = Mirror((0,1), 1)
mirror2 = Mirror((1,0.001), 1)
mirror3 = Mirror((-1,159), 1)
mirror4 = Mirror((150,1), 1)
table = TableOptique((20,20))
table.add(mirror3, (1.6,-5.5))
table.add(mirror1, (-5,4))
table.add(mirror2, (-8.7,1.69))
table.add(mirror4, (-3.81,-1.82))

table.draw()
table.draw_laser(rayon, (-1,1))
# print(mirror.reflect(rayon).k_vector)

plt.show()
