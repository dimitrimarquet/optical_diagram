#from msilib.schema import InstallUISequence
import numpy as np
import matplotlib.pyplot as plt
from test_wvl import wavelength_to_rgb

class Laser():
    """
    Représente un faisceau laser avec une longueur d'onde et un vecteur d'onde.

    Attributs:
        wavelength (float): Longueur d'onde du laser (en nm ou m selon l'unité choisie).
        k_vector (tuple): Vecteur d'onde (k_x, k_y) représentant la direction de propagation.
    """
    def __init__(self, wavelength, k_vector, stokes_vector = None): #, intensity = 1):
        self.wavelength = wavelength
        self.k_vector = k_vector
        if stokes_vector is None : 
            self.stokes = np.array([1,0,0,0]) 
        else :
            self.stokes = stokes_vector
    @property
    def intensity(self):
        return self.stokes[0] 
    
class Instruments():
    """
    Classe de base pour tous les instruments optiques.
    """
    def __init__(self, name):
        self.name = name
        pass

class Mirror(Instruments):
    """
    Représente un miroir optique avec une orientation et une réflectivité.

    Attributs:
        orientation (tuple): Vecteur d'orientation du miroir (mx, my).
        reflectivity (float): Coefficient de réflectivité du miroir (0 à 1).
    """

    def __init__(self, orientation_vector, reflectivity, name = 'Mirror'):
        super().__init__(name)
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
        
        #calculation of theta, angle of inclination of the mirror with respect to x
        if mx != 0:
            theta = -np.pi/2 + np.arctan(my/mx) 
        else :
            theta = 0

        #calculation of the reflected k vector k_out
        kx_theta = kx_init*np.cos(2*theta) + ky_init*np.sin(2*theta)
        ky_theta = - ky_init*np.cos(2*theta) + kx_init*np.sin(2*theta)
        k_out = (kx_theta, ky_theta)

        stokes_in = rayon.stokes.copy()

        stokes_out = stokes_in * self.reflectivity
        #print('mirror', stokes_in)
        #print('mirror', stokes_out)
        return Laser(wvl, k_out, stokes_vector=stokes_out)
   
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

        #calculation of theta, angle of inclination of the mirror with respect to x
        if mx != 0:
            theta = - np.pi/2 + np.arctan(my/mx)
        else :
            theta = 0
        
        #calculation of arrays (x_array, y_array) representing the mirror in (x,y) plane as a segment
        x_max = x_position + np.cos(theta)*length/2
        x_min = x_position - np.cos(theta)*length/2

        y_max = y_position + np.sin(theta)*length/2
        y_min = y_position - np.sin(theta)*length/2

        x_array = np.linspace(x_min, x_max, 10)
        y_array = np.linspace(y_min, y_max, 10)

        return [(x_array, y_array)]

class BeamSplitter(Instruments):
    """
    Représente un beamsplitter (lame séparatrice) qui divise le faisceau en deux.
    - Une partie est réfléchie (comme un miroir).
    - Une partie est transmise (traverse tout droit).
    
    Attributs:
        orientation (tuple): Vecteur d'orientation (mx, my).
        ratio (float): Ratio de transmission/réflexion de 0.5 pour un BS 50:50).
    """

    def __init__(self, orientation_vector, ratio=0.5):
        self.orientation = orientation_vector
        self.ratio = ratio  # 0.5 signifie 50% transmis, 50% réfléchis

    def reflect(self, rayon):
        """
        Gère la partie RÉFLÉCHIE du faisceau (Comportement identique au Miroir).
        """
        k_init = rayon.k_vector
        wvl = rayon.wavelength
        intensity = rayon.intensity
        mx, my = self.orientation
        kx_init, ky_init = k_init
        
        # Calcul de theta (angle d'inclinaison)
        if mx != 0:
            theta = -np.pi/2 + np.arctan(my/mx) 
        else :
            theta = 0

        # Calcul du vecteur k réfléchi (loi de la réflexion)
        kx_theta = kx_init*np.cos(2*theta) + ky_init*np.sin(2*theta)
        ky_theta = - ky_init*np.cos(2*theta) + kx_init*np.sin(2*theta)
        k_out = (kx_theta, ky_theta)

        coef_reflexion  = 1 - self.ratio 
        new_stokes = np.array(rayon.stokes.copy()) * coef_reflexion
        # On renvoie un NOUVEAU laser qui part dans la direction réfléchie
        return Laser(wvl, k_out, stokes_vector=new_stokes)

    def transmit(self, rayon):
        """
        Gère la partie TRANSMISE du faisceau.
        Dans un modèle simple (lame mince), le vecteur k ne change pas de direction.
        """

        new_stokes = np.array(rayon.stokes.copy())*self.ratio
        ###print("new_stokes", new_stokes)
        # Le faisceau continue tout droit : même vecteur k, même longueur d'onde
        return Laser(rayon.wavelength, rayon.k_vector, stokes_vector=new_stokes)
   
    def display(self, position, length=2):
        """
        Affiche un Cube Séparateur (BeamSplitter Cube).
        La 'length' définit la diagonale du cube (la surface active).
        Le cube est construit autour de cette diagonale.
        """
        xc, yc = position
        mx, my = self.orientation

        # 1. Calcul de l'angle (Theta)
        if mx != 0:
            theta = -np.pi/2 + np.arctan(my/mx)
        else:
            theta = 0

        # 2. Calcul des vecteurs
        # Rayon du cercle qui englobe le carré (demi-longueur de la diagonale)
        r = length / 2
        
        # Décalage pour la diagonale principale (Surface active)
        dx = r * np.cos(theta)
        dy = r * np.sin(theta)

        # 3. Calcul des 4 coins du carré
        # Coin 1 et 2 : Les extrémités du miroir (Diagonale active)
        p1 = (xc - dx, yc - dy)
        p2 = (xc + dx, yc + dy)
        
        # Coin 3 et 4 : Les autres coins (obtenus par rotation de 90° du vecteur dx,dy)
        # Vecteur perpendiculaire à (dx, dy) est (-dy, dx)
        p3 = (xc - dy, yc + dx)
        p4 = (xc + dy, yc - dx)

        # 4. Construction des tracés
        
        # SEGMENT A (Bleu) : La diagonale active (de P1 à P2)
        x_diag = np.array([p1[0], p2[0]])
        y_diag = np.array([p1[1], p2[1]])

        # SEGMENT B (Noir) : Le contour du carré
        # On relie les points dans l'ordre du périmètre : P1 -> P4 -> P2 -> P3 -> P1
        # (L'ordre P3/P4 dépend du sens trigo, mais tant qu'on fait le tour, c'est bon)
        x_box = np.array([p1[0], p4[0], p2[0], p3[0], p1[0]])
        y_box = np.array([p1[1], p4[1], p2[1], p3[1], p1[1]])

        return [(x_diag, y_diag), (x_box, y_box)]

class Polarizer(Instruments): 

    def __init__(self, orientation_vector , angle_deg):
        self.orientation = orientation_vector
        self.angle_deg = angle_deg 
        self.angle_rad = np.radians(angle_deg)

        theta = self.angle_rad
        c2 = np.cos(2*theta) 
        s2 = np.sin(2*theta)

        #Matrice de Mayer pour un polariseur linéaire 
        self.polarizer_matrix = 0.5 * np.array([
            [1,  c2,  s2,  0],
            [c2, c2**2, c2*s2, 0],
            [s2, c2*s2, s2**2, 0],
            [0,   0,   0,   0]       
        ])


    def transmit(self,laser):
        stokes_in = laser.stokes 
        stokes_out = np.dot(self.polarizer_matrix, stokes_in)
        # --- DEBUG : AFFICHER LE VECTEUR DE STOKES ---
        ###print(f"\n--- PASSAGE POLARISEUR ({self.angle_deg}°) ---")
        # np.round permet d'arrondir pour éviter les 0.000000001
        ###print(f"Stokes Entrée : {np.round(stokes_in, 3)}") 
        ###print(f"Stokes Sortie : {np.round(stokes_out, 3)}")
        ###print(f"Intensité transmise : {stokes_out[0]:.2f}")
        ###print("-----------------------------------------")
        # ---------------------------------------------
        return Laser(laser.wavelength, laser.k_vector, stokes_out)
    
    def reflect(self, laser):
        return None 
    
    def display(self,position, length = 2):
        x_position, y_position = position
        mx,my = self.orientation

        # dx = 0 
        # dy = length / 2
        
        # x_cadre = np.array([x_c, x_c])
        # y_cadre = np.array([y_c - dy, y_c + dy])
        
        # len_axis = length / 3
        # dx_a = (len_axis/2) * np.cos(self.angle_rad + np.pi/2) # +pi/2 car 0° = horizontal optique mais vertical sur le plot
        # dy_a = (len_axis/2) * np.sin(self.angle_rad + np.pi/2)
        
        # x_axis = np.array([x_c - dx_a, x_c + dx_a])
        # y_axis = np.array([y_c - dy_a, y_c + dy_a])
        
        # return [(x_axis, y_axis), (x_cadre, y_cadre)]

        #calculation of theta, angle of inclination of the mirror with respect to x
        if mx != 0:
            theta = - np.pi/2 + np.arctan(my/mx)
        else :
            theta = 0
        
        #calculation of arrays (x_array, y_array) representing the mirror in (x,y) plane as a segment
        x_max = x_position + np.cos(theta)*length/2
        x_min = x_position - np.cos(theta)*length/2

        y_max = y_position + np.sin(theta)*length/2
        y_min = y_position - np.sin(theta)*length/2

        x_array = np.linspace(x_min, x_max, 10)
        y_array = np.linspace(y_min, y_max, 10)

        return [(x_array, y_array), (x_array, y_array)]

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

    def remove(self, optics):
        if optics in self.keys():
            del self[optics]
        else :
            print('unknown optics')

    def draw(self, ax = None):
        """
        Affiche tous les instruments optiques présents sur la table.
        """
        if ax == None :
            for optics in self.keys():

                list_of_segments = optics.display(self[optics])

                if isinstance(optics, BeamSplitter): 
                    xd,yd = list_of_segments[0]
                    plt.plot(xd, yd, 'k-', linewidth=2, zorder = 10)

                    xc,yc = list_of_segments[1]
                    plt.plot(xc, yc, 'k-', linewidth = 2  )

                elif isinstance(optics, Mirror):
                    xm, ym = list_of_segments[0]
                    plt.plot(xm,ym, "k-", linewidth = 2, zorder = 10)

                elif isinstance(optics, Polarizer):
                    xc, yc = list_of_segments[0]
                    xa, ya = list_of_segments[1]

                    #  plt.plot(xc, yc, color="magenta", linewidth=3, alpha=0.5, label="Polarizer Body")
                    #  plt.plot(xa, ya, color="black", linewidth=2, label="Axis") 
                    plt.plot(xc, yc, color="magenta", linewidth=3, alpha=0.5)
                    plt.plot(xa, ya, color="black", linewidth=2) 
                custom_name = getattr(optics, 'name', None)

                if custom_name: 
                    label_text = custom_name
                else:
                    if isinstance(optics,Mirror): 
                        label_text = "Mirror"
                    elif isinstance(optics, BeamSplitter):
                        label_text = "BS"
                    elif isinstance(optics, Polarizer):
                        label_text = f"Pol {optics.angle_deg}°"
                print("draw_etx")
                x_cen, y_cen = self[optics]
                plt.text(x_cen + 1 ,y_cen + 1, label_text, color = 'black', fontsize = 9, ha = 'center', fontweight = 'bold')

        else:
            print('test')
            for optics in self.keys():
                list_of_segments = optics.display(self[optics])

                if isinstance(optics, BeamSplitter): 
                    xd,yd = list_of_segments[0]
                    ax.plot(xd, yd, 'k-', linewidth=2, zorder = 10)

                    xc,yc = list_of_segments[1]
                    ax.plot(xc, yc, 'k-', linewidth = 2  )

                elif isinstance(optics, Mirror):
                    xm, ym = list_of_segments[0]
                    ax.plot(xm,ym, "k-", linewidth = 2, zorder = 10)

                elif isinstance(optics, Polarizer):
                    xc, yc = list_of_segments[0]
                    xa, ya = list_of_segments[1]

                    #  plt.plot(xc, yc, color="magenta", linewidth=3, alpha=0.5, label="Polarizer Body")
                    #  plt.plot(xa, ya, color="black", linewidth=2, label="Axis") 
                    ax.plot(xc, yc, color="magenta", linewidth=3, alpha=0.5)
                    ax.plot(xa, ya, color="black", linewidth=2) 
            
                custom_name = getattr(optics, 'name', None)

                if custom_name: 
                    label_text = custom_name
                else:
                    if isinstance(optics,Mirror): 
                        label_text = "Mirror"
                    elif isinstance(optics, BeamSplitter):
                        label_text = "BS"
                    elif isinstance(optics, Polarizer):
                        label_text = f"Pol {optics.angle_deg}°"
                print("draw_etx")
                x_cen, y_cen = self[optics]
                plt.text(x_cen + 1 ,y_cen + 1, label_text, color = 'black', fontsize = 9, ha = 'center', fontweight = 'bold')

        #matplotlib adjustments
        xlim, ylim = self.size
        plt.xlim(-xlim/2, xlim/2)
        plt.ylim(-ylim/2, ylim/2)
        plt.gca().set_aspect('equal', adjustable='box')
        #plt.grid(True, linestyler = '--', alpha = 0.5)
    
    def path_laser(self, laser, position_init, segments=None, outputs = None):
        """
        Calcule le trajet d'un faisceau laser sur la table optique.

        Args:
            laser (Laser): Faisceau laser à propager.
            position_init (tuple): Position initiale (x, y) du faisceau.
            epsilon (float): Précision pour la détection des intersections.

        Returns:
            tuple: Deux listes (x_array, y_array) représentant le trajet du faisceau.
        """
        if segments is None:
            segments = []

        ray = laser
        x_init, y_init = position_init #position of the source of the laser
        x_size, y_size = self.size #dimensions of the table
        x_array = [x_init] #array of x coordinates needed for tracing the laser path
        y_array = [y_init] #array of y coordinates needed for tracing the laser path
        intensity_array = [ray.intensity]
        d_pos = np.sqrt(x_size**2 + y_size**2)*2 #max distance to calculate
        used_optics = [] #to avoid the program to loop on the same optics -> to modify for cavities !!

        i = 0
        running = True 

        while running and np.abs(x_array[-1]) < x_size/2 and np.abs(y_array[-1]) < y_size/2 and i < 100:

            #finding of the potential optics that could be touched by the laser considering its initial k vector, using the cone approach
            list_optics_in_cone = []
            kx_init, ky_init = ray.k_vector
            norm_init = np.sqrt(kx_init**2 + ky_init**2)
            
            if ray.stokes[0] < 0.01:
                return segments

            ###print('INTENSITY', i, ray.intensity)
            for optics in self.keys():

                if len(used_optics) > 0 and optics == used_optics[-1]:
                    continue

                x_optics, y_optics = self[optics]
                segments_display = optics.display(self[optics])
                x_array_optics, y_array_optics =segments_display[0]

                i_xmax_optics = np.argmax(x_array_optics)
                i_xmin_optics = np.argmin(x_array_optics)
                i_ymax_optics = np.argmax(y_array_optics)
                i_ymin_optics = np.argmin(y_array_optics)

                if i_xmax_optics != i_xmin_optics:
                    kmax_x, kmax_y = (x_array_optics[i_xmax_optics] - x_init, y_array_optics[i_xmax_optics] - y_init)
                    kmin_x, kmin_y = (x_array_optics[i_xmin_optics] - x_init, y_array_optics[i_xmin_optics] - y_init)
                else :
                    kmax_x, kmax_y = (x_array_optics[i_ymax_optics] - x_init, y_array_optics[i_ymax_optics] - y_init)
                    kmin_x, kmin_y = (x_array_optics[i_ymin_optics] - x_init, y_array_optics[i_ymin_optics] - y_init)


                norm_max = np.sqrt(kmax_x**2 + kmax_y**2)
                norm_min = np.sqrt(kmin_x**2 + kmin_y**2)

                #test to be on the right side of the cone if scalar product positive
                scalar_max = np.dot(np.array([kmax_x, kmax_y])/norm_max, np.array([kx_init, ky_init])/norm_init)
                scalar_min = np.dot(np.array([kmin_x, kmin_y])/norm_min, np.array([kx_init, ky_init])/norm_init)

                #test to be in the cone if the angle with respect to x is between the two
                test_max = np.arccos(min(scalar_max, 1)) < np.arccos(min(np.dot(np.array([kmax_x, kmax_y])/norm_max, np.array([kmin_x, kmin_y])/norm_min), 1))
                test_min = np.arccos(min(scalar_min, 1)) < np.arccos(min(np.dot(np.array([kmax_x, kmax_y])/norm_max, np.array([kmin_x, kmin_y])/norm_min), 1))
                print('test', np.array([kmax_x, kmax_y])/norm_max)
                print('test', np.array([kmin_x, kmin_y])/norm_min)

                test_optics_in_cone = scalar_max > 0 and scalar_min > 0 and test_max and test_min

                if test_optics_in_cone and optics not in used_optics[-1:]:
                    list_optics_in_cone.append(optics)
                    print(i, list_optics_in_cone)
                    print(ky_init, kx_init)

            
            if len(list_optics_in_cone) != 0:

                #discrimination of the optics that is actually going to be touched (closest from the initial position)
                list_position_optics_in_cone = [self[el] for el in list_optics_in_cone]
                list_distance_optics_in_cone = [np.sqrt((el[0] - x_array[-1])**2 + (el[1] - y_array[-1])**2) for el in list_position_optics_in_cone]

                optics = list_optics_in_cone[np.argmin(list_distance_optics_in_cone)]
                position_optics = self[optics]

                #calculation of the intersection point and reflection except in the case where the laser comes on the laser side
                x_optics, y_optics = position_optics
                segments_display = optics.display(position_optics)
                x_array_optics, y_array_optics = segments_display[0]

                i_xmax_optics = np.argmax(x_array_optics)
                i_xmin_optics = np.argmin(x_array_optics)
                i_ymax_optics = np.argmax(y_array_optics)
                i_ymin_optics = np.argmin(y_array_optics)

                if i_xmax_optics != i_xmin_optics:
                    kmax_x, kmax_y = (x_array_optics[i_xmax_optics] - x_init, y_array_optics[i_xmax_optics] - y_init)
                    kmin_x, kmin_y = (x_array_optics[i_xmin_optics] - x_init, y_array_optics[i_xmin_optics] - y_init)
                else :
                    kmax_x, kmax_y = (x_array_optics[i_ymax_optics] - x_init, y_array_optics[i_ymax_optics] - y_init)
                    kmin_x, kmin_y = (x_array_optics[i_ymin_optics] - x_init, y_array_optics[i_ymin_optics] - y_init)

                norm_max = np.sqrt(kmax_x**2 + kmax_y**2)
                norm_min = np.sqrt(kmin_x**2 + kmin_y**2)

                if x_array_optics[i_xmax_optics] - x_array_optics[i_xmin_optics] != 0: #if the mirror is not along y
                    coef_miroir = (y_array_optics[i_xmax_optics] - y_array_optics[i_xmin_optics])/(x_array_optics[i_xmax_optics] - x_array_optics[i_xmin_optics])
                    coef_laser = ky_init/kx_init

                    if coef_miroir-coef_laser != 0: #if the laser does not land on the mirror's side
                        
                        x_intersect = (y_init-y_optics + coef_miroir*x_optics - coef_laser*x_init)/(coef_miroir-coef_laser)
                        y_intersect = coef_miroir*(x_intersect-x_optics) + y_optics
                        x_array.append(x_intersect)
                        y_array.append(y_intersect)

                        if isinstance(optics, Mirror): 
                            ###print('testjpp', ray.intensity)
                            x_init, y_init = x_intersect, y_intersect
                            ray = optics.reflect(ray)

                        #Prise en compte du Beam Splitter
                        elif isinstance(optics, BeamSplitter): 
                            transmitted_ray = optics.transmit(ray) #calcul du rayon transmis

                            current_intensity_BS = transmitted_ray.intensity
                            segments.append((list(x_array), list(y_array), current_intensity_BS)) #sauvergard du chemin actuel avant de partir 
                            self.path_laser(transmitted_ray, (x_intersect,y_intersect),segments, outputs = outputs) #Relance de la fonction pour le rayon transmis
                            x_init, y_init = x_intersect, y_intersect
                            ray = optics.reflect(ray)                            

                        elif isinstance(optics, Polarizer):
                            ##print('test_pol')
                            current_intensity_pol = ray.intensity
                            segments.append((list(x_array),list(y_array),current_intensity_pol))
                            ray = optics.transmit(ray)

                            if ray.stokes[0] < 0.01:                                
                                return segments
                            # else : 
                            #     epsilon = 0.1 
                            #     kx, ky = ray.k_vector 
                            #     norm_k = np.sqrt(kx**2 + ky**2)
                            #     # x_init = x_intersect + (kx/norm_k) * epsilon
                            #     # y_init = y_intersect + (ky/norm_k) * epsilon

                            #     # x_array = [x_init]
                            #     # y_array = [y_init]

                            #     used_optics.append(optics)
                        
                        intensity_array.append(ray.intensity)
                        used_optics.append(optics)
                    
                    else :
                        x_array.append(x_optics)
                        y_array.append(y_optics)
                        intensity_array.append(ray.intensity)
                        return (x_array, y_array)

                else : #if the mirror is along y
                    used_optics.append(optics)
                    ###print("not implemented yet...")
                        
            else:
                x_array.append(x_array[-1] + d_pos*kx_init)
                y_array.append(y_array[-1] + d_pos*ky_init)
                intensity_array.append(ray.intensity)

                if outputs is not None: 
                    outputs.append(ray)
                    print('test output')
                # running = False

            i += 1
            
        segments.append((x_array, y_array, ray.intensity))  
        return segments
    
    def draw_laser(self, laser, position_source, ax = None):
        """
        Affiche le trajet d'un faisceau laser sur la table optique.

        Args:
            laser (Laser): Faisceau laser à afficher.
            position_source (tuple): Position initiale (x, y) du faisceau.
        """
        
        all_segments = self.path_laser(laser, position_source)
        
        color = wavelength_to_rgb(laser.wavelength)
        r = color[0]
        g = color[1]
        b = color[2]
        k = 0
        for i, (x_array, y_array, intensity) in enumerate(all_segments):
            valeur_intensity = float(intensity)
            if valeur_intensity < 0.01:  
                continue
            alpha_val = min(intensity, 1.0)
            if ax == None:
                plt.plot(x_array, y_array, color = (r/256,g/256,b/256), linewidth = 1)
            else :
                ax.plot(x_array, y_array, color = (r/256,g/256,b/256), linewidth = 1)
            k+=1

    def report(self, laser,position_source, filename ="rapport_optique.txt"):
        """
        Génère un fichier texte contenant les caractéristiques de tous les faisceaux
        qui sortent de la table optique.

        """
        final_beam = []

        self.path_laser(laser, position_source,segments=[], outputs = final_beam)

        with open(filename, 'w', encoding="utf-8") as f:
            f.write("===RAPPORT DE SIMULATION OPTIQUE ===\n")
            f.write(f"Source Lambda : {laser.wavelength} nm\n")
            f.write(f"Source Stokes : {laser.stokes}\n")
            f.write("-" * 40 + "\n\n")
            print(len(final_beam))
            if not final_beam:
                f.write("Aucune faisceau n'est sorti de la table")

            for i, ray in enumerate(final_beam):
                f.write(f"--- FAISCEAU DE SORTIE #{i+1} ---\n")

                kx, ky = ray.k_vector 
                angle_rad = np.arctan2(ky,kx)
                angle_def = np.degrees(angle_rad)

                I, Q, U, V = ray.stokes 
                pol_type = "Inconnue"

                polarized_intensity = np.sqrt(Q**2 + U**2 + V**2)

                if I > 1e-9:
                    dop = polarized_intensity / I
                else : 
                    dop = 0 
                
                if dop < 0.1: 
                    pol_type = "Non polarisée"
                elif abs(V) > 0.1 and abs(Q) < 0.1 and abs(U) < 0.1: 
                    pol_type = "Circulaire"
                elif abs(V) < 0.1:
                    pol_type = "Linéaire"
                else : 
                    pol_type = "Elliptique"

                f.write(f"Intensité Finale : {I:.4f} ({(I/laser.intensity)*100:.1f}% de l'entrée)\n")
                f.write(f"Direction (k)    : ({kx:.3f}, {ky:.3f})\n")
                f.write(f"Angle            : {angle_def:.2f} deg\n")
                f.write(f"Vecteur Stokes   : [{I:.3f}, {Q:.3f}, {U:.3f}, {V:.3f}]\n")
                f.write(f"Type Polarisation: {pol_type}\n")
                f.write("\n")
                    
        print(f"Rapport généré : {filename}")




# k = (1,0)
# # k = (2.1,1)
# # k2 = (-2.1, -1)

# rayon = Laser(633, k, stokes_vector=[1,1,0,0])
# # rayon2 = Laser(600, k2)
# # rayon = Laser(620, k)
# mirror1 = Mirror((0,1), 1)
# # mirror2 = Mirror((1,0.001), 1)
# # mirror3 = Mirror((-1,159), 1)
# # mirror4 = Mirror((150,1), 1)

# table = TableOptique((20,20))


# table.add(mirror1, (-5,4))
# # table.add(mirror2, (-8.7,1.69))
# # table.add(mirror4, (-3.81,-1.82))

# #test_cavity

# mirror1 = Mirror((1,2), 1)
# mirror2 = Mirror((5,1), 1)
# mirror3 = Mirror((1,1), 1)
# mirror4 = Mirror((1,10), 1)
# mirror5 = Mirror((1,0.01), 1)
# mirror6 = Mirror((-0.1,1), 1)

# BS1 = BeamSplitter((1,1.1),ratio=0.3)
# P2 = Polarizer((0.1,1),90) 
# P4 = Polarizer((1,0.01),-45) 

# P3 = Mirror((0.25,1),1) 

# # laser_out = P2.transmit(rayon)
# # ###print(f"Intensité après Polariseur 90° : {laser_out.intensity:.5f}")

# table.add(mirror1, (7,0))
# table.add(mirror2, (1.8,-6.5))
# table.add(mirror4, (8.55,-2.2))
# table.add(mirror5, (5,0))
# table.add(mirror6, (1,6))

# table.add(BS1,(1,0.1))
# table.add(P2,(1.5,-5))
# table.add(P3,(1.7,-8))
# table.add(P4,(5.79, -4.79))

# # mirror1 = Mirror((1,2), 1)
# # mirror3 = Mirror((1,2), 1)

# # mirror2 = Mirror((1,-1), 1)
# # mirror4 = Mirror((50,0.1), 1)
# # table.add(mirror1, (7,0))
# # table.add(mirror3, (7.5,0.5))
# # table.add(mirror2, (4.82,3))
# # table.add(mirror4, (8.93,-2.58))


# table.draw()
# table.draw_laser(rayon, (0,0))



# plt.show()