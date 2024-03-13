from random import random, uniform, seed
from matplotlib import pyplot as plt
import numpy as np
import math
import os
#from tqdm import trange

from geom import distance, intersection, parallel, angle, orthogonal, solve
from geom import Point, Vector, Line, Plane

from funcs import generate_n_planes, get_distribution_form_r_k
from funcs import Cube, RegularTriangularPrism, Parallelepiped
from funcs import make_random_plane_into_cube, make_random_plane_into_prism, make_random_plane_into_tetra

from funcs import WCCoPrism, make_random_plane_into_wcco

from funcs import area_of_polygon_3d, get_angles, secant_polygon, characteristic_of_section, norm_dict
from funcs import rotate_points_to_Oxy, is_convex_polygon_2D, perimeter_3D, semiaxes_from_polygon_3D, minmax_feret

#from funcs import make_distributions

from PyQt6.QtCore import QDateTime, Qt, QTimer, QThread, pyqtSignal
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDateTimeEdit,
        QDial, QDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QSpinBox, QStyleFactory, QTableWidget, QTabWidget, QTextEdit,
        QVBoxLayout, QWidget, QFormLayout, QDoubleSpinBox, QFileDialog)

from PyQt6.QtGui import QPainter, QPalette

import random

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class WorkerThread(QThread):
    update_signal = pyqtSignal(int)
    

    def make_distributions(self, polyhedron, way, PATH, type_figure='parallelepiped', n=10000, params=None):
        '''
        генерирует n плоскостей способом way и возвращает плотность распределения
        геометрических характеристик сечений: углов, площадей, 
        '''
        #углы
        angles = [0]*181

        #ПАРАЛЛЕЛЕПИПЕД
        if type_figure == 'parallelepiped':

            if params:
                a_side = params['side_a']
                b_side = params['side_b']
                c_side = params['side_c']
            else:
                a_side = 1.
                b_side = 2.
                c_side = 3.
                params = {'side_a': 1.,
                          'side_b': 2.,
                          'side_c': 3.}

            #path_ = PATH + r'\parallelepiped_a=' + str(a_side) + 'b=' + str(b_side) + 'c=' + str(c_side)
            path_ = os.path.join(PATH, f"parallelepiped_a={a_side}b={b_side}c={c_side}")

            try:
                angles = np.loadtxt(os.path.join(path_, f'anlges_N={n}.txt'))
                areas = np.load(os.path.join(path_, f'areas_N={n}.npy'), allow_pickle=True).item()
                a_semi = np.load(os.path.join(path_, f'asemi_N={n}.npy'), allow_pickle=True).item()
                b_semi = np.load(os.path.join(path_, f'bsemi_N={n}.npy'), allow_pickle=True).item()
                perimeter = np.load(os.path.join(path_, f'perimeter_N={n}.npy'), allow_pickle=True).item()
                minferet = np.load(os.path.join(path_, f'minferet_N={n}.npy'), allow_pickle=True).item()
                maxferet = np.load(os.path.join(path_, f'maxferet_N={n}.npy'), allow_pickle=True).item()
                aspect = np.load(os.path.join(path_, f'aspect_N={n}.npy'), allow_pickle=True).item()
                sphericity = np.load(os.path.join(path_, f'sphericity_N={n}.npy'), allow_pickle=True).item()

            except:
                areas_dict = {} #пустой словать с ключами от 0.00 до a^2+b^2+c^2 #можно меньше
                step = 0.01
                current_key = 0.00
                while current_key <= 10*(a_side**2 + b_side**2 + c_side**2):
                    areas_dict.setdefault(round(current_key, 2), 0)
                    current_key += step


                perim_dict = {} #пустой словать с ключами от 0.00 до 6*max(a,b,c)
                step = 0.01
                current_key = 0.00
                while current_key <= 20*max(a_side, b_side, c_side):
                    perim_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                #max полуось <= max(a,b,c)
                #semiaxes = np.zeros((10*round(max(a_side, b_side, c_side)*100)+1, 10*round(max(a_side, b_side, c_side)*100)+1))

                a_semi_dict = {} #пустой словать с ключами от 0.00 до max(a,b,c) с шагом 0.01
                b_semi_dict = {}
                step = 0.01
                current_key = 0.00
                while current_key <= 10*max(a_side, b_side, c_side):
                    a_semi_dict.setdefault(round(current_key, 2), 0)
                    b_semi_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                minferet_dict = {} #пустой словать с ключами от 0.00 до max(a,b,c) с шагом 0.01
                maxferet_dict = {}
                aspect_dict = {}
                sphericity_dict = {}

                step = 0.01
                current_key = 0.00
                while current_key <= 3*a_side*b_side*c_side:
                    minferet_dict.setdefault(round(current_key, 2), 0)
                    maxferet_dict.setdefault(round(current_key, 2), 0)
                    current_key += step
                    
                current_key = 0.00    
                while current_key <= 1.:
                    aspect_dict.setdefault(round(current_key, 2), 0)
                    sphericity_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                num_errors = 0   


                for i in range(n):
                    sec_polygon = secant_polygon(polyhedron, way, params=params)
                    if len(sec_polygon) > 2:
                        polygonOxy = rotate_points_to_Oxy(sec_polygon)[0]
                        #check convexity
                        if is_convex_polygon_2D(polygonOxy):
                            ang_sec = get_angles(sec_polygon)
                            area = area_of_polygon_3d(sec_polygon)
                            perim = perimeter_3D(sec_polygon)
                            a_semi, b_semi = semiaxes_from_polygon_3D(sec_polygon)
                            
                            min_feret, max_feret = minmax_feret(polygonOxy)
                            aspect_ratio = min_feret/max_feret
                            sphericity = 2*(np.pi*area)**0.5/perim
                            if aspect_ratio <= 1.:
                                for ang in ang_sec:
                                    angles[ang] += 1
                                areas_dict[round(area, 2)] += 1
                                perim_dict[round(perim, 2)] += 1
                                a_semi_dict[round(a_semi, 2)] += 1
                                b_semi_dict[round(b_semi, 2)] += 1
                                
                                minferet_dict[round(min_feret, 2)] += 1
                                maxferet_dict[round(max_feret, 2)] += 1
                                aspect_dict[round(aspect_ratio, 2)] += 1
                                sphericity_dict[round(sphericity, 2)] += 1
                            else:
                                num_errors += 1
                        else:
                            num_errors += 1
                    else:
                        num_errors += 1

                    if i % (n//100) == 0:
                        self.update_signal.emit(int((i + 1)/n*100))
                
                print('precision is ', 1-num_errors/n)

                angles = np.array(list(np.array(angles)/np.sum(angles)))
                areas = norm_dict(areas_dict, step=0.01)
                perimeter = norm_dict(perim_dict, step=0.01)
                a_semi = norm_dict(a_semi_dict, step=0.01)
                b_semi = norm_dict(b_semi_dict, step=0.01)

                minferet = norm_dict(minferet_dict, step=0.01)
                maxferet = norm_dict(maxferet_dict, step=0.01)
                aspect = norm_dict(aspect_dict, step=0.01)
                sphericity = norm_dict(sphericity_dict, step=0.01)

                try:
                    os.mkdir(path_)
                except:
                    pass

                np.savetxt(os.path.join(path_, f'anlges_N={n}.txt'), angles)
                np.save(os.path.join(path_, f'areas_N={n}.npy'), areas)
                np.save(os.path.join(path_, f'perimeter_N={n}.npy'), perimeter)
                np.save(os.path.join(path_, f'asemi_N={n}.npy'), a_semi)
                np.save(os.path.join(path_, f'bsemi_N={n}.npy'), b_semi)
                #np.save(os.path.join(path_, f'semiaxes_N={n}.npy'), semiaxes)

                np.save(os.path.join(path_, f'minferet_N={n}.npy'), minferet)
                np.save(os.path.join(path_, f'maxferet_N={n}.npy'), maxferet)
                np.save(os.path.join(path_, f'aspect_N={n}.npy'), aspect)
                np.save(os.path.join(path_, f'sphericity_N={n}.npy'), sphericity)



        if type_figure == 'triangular prism':
            #[0.00, 0.01, 0.02, ..., 1.40, 1.41, 1.42]#
            #[0, 1, 2, ..., 140, 141, 142]

            if params:
                side = params['side']
            else:
                side = 1.
                params = {'side': 1.}

            path_ = os.path.join(PATH, f"tri-prism_side={side}")

            try:
                angles = np.loadtxt(os.path.join(path_, f'anlges_N={n}.txt'))
                areas = np.load(os.path.join(path_, f'areas_N={n}.npy'), allow_pickle=True).item()
                a_semi = np.load(os.path.join(path_, f'asemi_N={n}.npy'), allow_pickle=True).item()
                b_semi = np.load(os.path.join(path_, f'bsemi_N={n}.npy'), allow_pickle=True).item()
                perimeter = np.load(os.path.join(path_, f'perimeter_N={n}.npy'), allow_pickle=True).item()
                minferet = np.load(os.path.join(path_, f'minferet_N={n}.npy'), allow_pickle=True).item()
                maxferet = np.load(os.path.join(path_, f'maxferet_N={n}.npy'), allow_pickle=True).item()
                aspect = np.load(os.path.join(path_, f'aspect_N={n}.npy'), allow_pickle=True).item()
                sphericity = np.load(os.path.join(path_, f'sphericity_N={n}.npy'), allow_pickle=True).item()

            except:
                areas = [0]*143

                areas_dict = {} #пустой словать с ключами от 0.00 до 1.42 с шагом 0.01
                step = 0.01
                current_key = 0.00
                while current_key <= 3**0.5*side*side:
                    areas_dict.setdefault(round(current_key, 2), 0)
                    current_key += step


                perim_dict = {} #пустой словать с ключами от 0.00 до 10.00 с шагом 0.01
                step = 0.01
                current_key = 0.00
                while current_key <= 4*side*2:
                    perim_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                #semiaxes = np.zeros((round(2*side*100), round(2*side*100)), dtype=int)

                a_semi_dict = {} #пустой словать с ключами от 0.00 до 3.00 с шагом 0.01
                b_semi_dict = {}
                step = 0.01
                current_key = 0.00
                while current_key <= 2*side:
                    a_semi_dict.setdefault(round(current_key, 2), 0)
                    b_semi_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                minferet_dict = {} #пустой словать с ключами от 0.00 до max(a,b,c) с шагом 0.01
                maxferet_dict = {}
                aspect_dict = {}
                sphericity_dict = {}

                step = 0.01
                current_key = 0.00
                while current_key <= 5*side:
                    minferet_dict.setdefault(round(current_key, 2), 0)
                    maxferet_dict.setdefault(round(current_key, 2), 0)
                    current_key += step
                    
                current_key = 0.00    
                while current_key <= 1.:
                    aspect_dict.setdefault(round(current_key, 2), 0)
                    sphericity_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                num_errors = 0   


                for i in range(n):
                    sec_polygon = secant_polygon(polyhedron, way, params=params)
                    if len(sec_polygon) > 2:
                        polygonOxy = rotate_points_to_Oxy(sec_polygon)[0]
                        #check convexity
                        if is_convex_polygon_2D(polygonOxy):
                            ang_sec = get_angles(sec_polygon)
                            area = area_of_polygon_3d(sec_polygon)
                            perim = perimeter_3D(sec_polygon)
                            a_semi, b_semi = semiaxes_from_polygon_3D(sec_polygon)
                            
                            min_feret, max_feret = minmax_feret(polygonOxy)
                            aspect_ratio = min_feret/max_feret
                            sphericity = 2*(np.pi*area)**0.5/perim
                            if aspect_ratio <= 1.:
                                for ang in ang_sec:
                                    angles[ang] += 1
                                areas_dict[round(area, 2)] += 1
                                perim_dict[round(perim, 2)] += 1
                                a_semi_dict[round(a_semi, 2)] += 1
                                b_semi_dict[round(b_semi, 2)] += 1
                                
                                minferet_dict[round(min_feret, 2)] += 1
                                maxferet_dict[round(max_feret, 2)] += 1
                                aspect_dict[round(aspect_ratio, 2)] += 1
                                sphericity_dict[round(sphericity, 2)] += 1
                            else:
                                num_errors += 1
                        else:
                            num_errors += 1
                    else:
                        num_errors += 1
                    if i % (n//100) == 0:
                        self.update_signal.emit(int((i + 1)/n*100))


                print('precision is ', 1-num_errors/n)

                angles = np.array(list(np.array(angles)/np.sum(angles)))
                areas = norm_dict(areas_dict, step=0.01)
                perimeter = norm_dict(perim_dict, step=0.01)
                a_semi = norm_dict(a_semi_dict, step=0.01)
                b_semi = norm_dict(b_semi_dict, step=0.01)

                minferet = norm_dict(minferet_dict, step=0.01)
                maxferet = norm_dict(maxferet_dict, step=0.01)
                aspect = norm_dict(aspect_dict, step=0.01)
                sphericity = norm_dict(sphericity_dict, step=0.01)

                try:
                    os.mkdir(path_)
                except:
                    pass

                np.savetxt(os.path.join(path_, f'anlges_N={n}.txt'), angles)
                np.save(os.path.join(path_, f'areas_N={n}.npy'), areas)
                np.save(os.path.join(path_, f'perimeter_N={n}.npy'), perimeter)
                np.save(os.path.join(path_, f'asemi_N={n}.npy'), a_semi)
                np.save(os.path.join(path_, f'bsemi_N={n}.npy'), b_semi)
                #np.save(os.path.join(path_, f'semiaxes_N={n}.npy'), semiaxes)

                np.save(os.path.join(path_, f'minferet_N={n}.npy'), minferet)
                np.save(os.path.join(path_, f'maxferet_N={n}.npy'), maxferet)
                np.save(os.path.join(path_, f'aspect_N={n}.npy'), aspect)
                np.save(os.path.join(path_, f'sphericity_N={n}.npy'), sphericity)



        if type_figure == 'hex prism':

            r = params['r']
            k = params['k']

            path_ = os.path.join(PATH, f"hex-prism_r={r}k={k}")

            try:
                angles = np.loadtxt(os.path.join(path_, f'anlges_N={n}.txt'))
                areas = np.load(os.path.join(path_, f'areas_N={n}.npy'), allow_pickle=True).item()
                a_semi = np.load(os.path.join(path_, f'asemi_N={n}.npy'), allow_pickle=True).item()
                b_semi = np.load(os.path.join(path_, f'bsemi_N={n}.npy'), allow_pickle=True).item()
                perimeter = np.load(os.path.join(path_, f'perimeter_N={n}.npy'), allow_pickle=True).item()
                minferet = np.load(os.path.join(path_, f'minferet_N={n}.npy'), allow_pickle=True).item()
                maxferet = np.load(os.path.join(path_, f'maxferet_N={n}.npy'), allow_pickle=True).item()
                aspect = np.load(os.path.join(path_, f'aspect_N={n}.npy'), allow_pickle=True).item()
                sphericity = np.load(os.path.join(path_, f'sphericity_N={n}.npy'), allow_pickle=True).item()

            except:
                angles = [0]*181

                areas_dict = {} #пустой словать с ключами от 0.00 до 1.42 с шагом 0.01
                step = 0.01
                current_key = 0.00
                while current_key <= 10.:
                    areas_dict.setdefault(round(current_key, 2), 0)
                    current_key += step


                perim_dict = {} #пустой словать с ключами от 0.00 до 10.00 с шагом 0.01
                step = 0.01
                current_key = 0.00
                while current_key <= 10.:
                    perim_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                #semiaxes = np.zeros((round(k*2*60)+20, round(k*2*60)+20), dtype=int)

                a_semi_dict = {} #пустой словать с ключами от 0.00 до 3.00 с шагом 0.01
                b_semi_dict = {}
                step = 0.01
                current_key = 0.00
                while current_key <= round(k*2*60)+20:
                    a_semi_dict.setdefault(round(current_key, 2), 0)
                    b_semi_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                minferet_dict = {} #пустой словать с ключами от 0.00 до max(a,b,c) с шагом 0.01
                maxferet_dict = {}
                aspect_dict = {}
                sphericity_dict = {}

                step = 0.01
                current_key = 0.00
                while current_key <= round(k*2*60) + 20.:
                    minferet_dict.setdefault(round(current_key, 2), 0)
                    maxferet_dict.setdefault(round(current_key, 2), 0)
                    current_key += step
                    
                current_key = 0.00    
                while current_key <= 1.:
                    aspect_dict.setdefault(round(current_key, 2), 0)
                    sphericity_dict.setdefault(round(current_key, 2), 0)
                    current_key += step

                num_errors = 0   


                for i in range(n):
                    sec_polygon = secant_polygon(polyhedron, way, params=params)
                    if len(sec_polygon) > 2:
                        polygonOxy = rotate_points_to_Oxy(sec_polygon)[0]
                        #check convexity
                        if is_convex_polygon_2D(polygonOxy):
                            ang_sec = get_angles(sec_polygon)
                            area = area_of_polygon_3d(sec_polygon)
                            perim = perimeter_3D(sec_polygon)
                            a_semi, b_semi = semiaxes_from_polygon_3D(sec_polygon)
                            
                            min_feret, max_feret = minmax_feret(polygonOxy)
                            aspect_ratio = min_feret/max_feret
                            sphericity = 2*(np.pi*area)**0.5/perim
                            if aspect_ratio <= 1.:
                                for ang in ang_sec:
                                    angles[ang] += 1
                                areas_dict[round(area, 2)] += 1
                                perim_dict[round(perim, 2)] += 1
                                a_semi_dict[round(a_semi, 2)] += 1
                                b_semi_dict[round(b_semi, 2)] += 1
                                
                                minferet_dict[round(min_feret, 2)] += 1
                                maxferet_dict[round(max_feret, 2)] += 1
                                aspect_dict[round(aspect_ratio, 2)] += 1
                                sphericity_dict[round(sphericity, 2)] += 1
                            else:
                                num_errors += 1
                        else:
                            num_errors += 1
                    else:
                        num_errors += 1
                    if i % (n//100) == 0:
                        self.update_signal.emit(int((i + 1)/n*100))


                print('precision is ', 1-num_errors/n)

                angles = np.array(list(np.array(angles)/np.sum(angles)))
                areas = norm_dict(areas_dict, step=0.01)
                perimeter = norm_dict(perim_dict, step=0.01)
                a_semi = norm_dict(a_semi_dict, step=0.01)
                b_semi = norm_dict(b_semi_dict, step=0.01)

                minferet = norm_dict(minferet_dict, step=0.01)
                maxferet = norm_dict(maxferet_dict, step=0.01)
                aspect = norm_dict(aspect_dict, step=0.01)
                sphericity = norm_dict(sphericity_dict, step=0.01)

                try:
                    os.mkdir(path_)
                except:
                    pass

                np.savetxt(os.path.join(path_, f'anlges_N={n}.txt'), angles)
                np.save(os.path.join(path_, f'areas_N={n}.npy'), areas)
                np.save(os.path.join(path_, f'perimeter_N={n}.npy'), perimeter)
                np.save(os.path.join(path_, f'asemi_N={n}.npy'), a_semi)
                np.save(os.path.join(path_, f'bsemi_N={n}.npy'), b_semi)
                #np.save(os.path.join(path_, f'semiaxes_N={n}.npy'), semiaxes)

                np.save(os.path.join(path_, f'minferet_N={n}.npy'), minferet)
                np.save(os.path.join(path_, f'maxferet_N={n}.npy'), maxferet)
                np.save(os.path.join(path_, f'aspect_N={n}.npy'), aspect)
                np.save(os.path.join(path_, f'sphericity_N={n}.npy'), sphericity)


        return angles, areas, perimeter, a_semi, b_semi, minferet, maxferet, aspect, sphericity



class CustomCanvas(FigureCanvas):
    def __init__(self, parent=None, size=5, dpi=100):
        self.figure = Figure(figsize=(size, size), dpi=dpi)
        self.axes = self.figure.add_subplot(111)
        super().__init__(self.figure)
        self.setParent(parent)
        
    
    def plot_distribution(self, type_figure, characteristic, PATH, N, thread, params=None):
        
        self.axes.clear()
        
        #просто нужна эта функция
        def extract_and_trim_dict(dictionary):
            keys = []
            values = []
            for key, value in dictionary.items():
                keys.append(key)
                values.append(value)

            # Обрезаем конец списков, пока последний элемент не равен 0.0
            while values and values[-1] == 0.0:
                del keys[-1]
                del values[-1]
            return keys, values
        
        
        
        if type_figure =='parallelepiped':
            if not params:
                params = {'side_a': 1.,
                          'side_b': 2.,
                          'side_c': 3.}
            type_figure += ' a=' + str(params['side_a']) + ' b=' + str(params['side_b']) + ' c=' + str(params['side_c'])
            
            angles, areas, perimeter, a_semi, b_semi, minferet, maxferet, aspect, sphericity = thread.make_distributions(Parallelepiped(**params), 
                                                                                       make_random_plane_into_tetra,
                                                                                       PATH=PATH,
                                                                                       type_figure='parallelepiped', 
                                                                                       n=N,
                                                                                       params=params)
            
        #ТРЕУГОЛЬНАЯ ПРИЗМА                
        if type_figure =='triangular prism':
            if not params:
                params = {'side': 1.}
            type_figure += ' height=' + str(params['side'])
            
            
            angles, areas, perimeter, a_semi, b_semi, minferet, maxferet, aspect, sphericity = thread.make_distributions(RegularTriangularPrism(**params), 
                                                                                    make_random_plane_into_prism, 
                                                                                    PATH=PATH,
                                                                                    type_figure='triangular prism', 
                                                                                    n=N,
                                                                                    params=params)


        if type_figure == 'hex prism':
            if not params:
                params = {'r': 1.,
                          'k': 1.}
            type_figure += ' r=' + str(params['r']) + ' k=' + str(params['k'])
            
            angles, areas, perimeter, a_semi, b_semi, minferet, maxferet, aspect, sphericity = thread.make_distributions(WCCoPrism(**params), 
                                                                                    make_random_plane_into_wcco, 
                                                                                    PATH=PATH,
                                                                                    type_figure='hex prism', 
                                                                                    n=N,
                                                                                    params=params)
            
            
        if characteristic =='angle':
            self.axes.plot(range(181), angles)
            self.axes.set_xlabel('angle')
            self.axes.set_ylabel('p')
            self.axes.set_title('Angles of ' + type_figure)
            self.draw()
            
        if characteristic =='area':
            self.axes.plot(*extract_and_trim_dict(areas))
            self.axes.set_xlabel('area')
            self.axes.set_ylabel('p')
            self.axes.set_title('Area of ' + type_figure)
            self.draw()
        
        if characteristic =='perimeter':
            self.axes.plot(*extract_and_trim_dict(perimeter))
            self.axes.set_xlabel('perimeter')
            self.axes.set_ylabel('p')
            self.axes.set_title('Perimeter of ' + type_figure)
            self.draw()

        if characteristic =='semiaxes':
            self.axes.plot(*extract_and_trim_dict(a_semi), label='semi-minor axis')
            self.axes.plot(*extract_and_trim_dict(b_semi), label='semi-major axis')
            self.axes.set_xlabel('semiaxes')
            self.axes.set_ylabel('p')
            self.axes.set_title('Semiaxes of ' + type_figure)
            self.axes.legend()
            self.draw()
        
        if characteristic =='Feret diameter':
            self.axes.plot(*extract_and_trim_dict(minferet), label='min Feret diameter')
            self.axes.plot(*extract_and_trim_dict(maxferet), label='max Feret diameter')
            self.axes.set_xlabel('Feret diameter')
            self.axes.set_ylabel('p')
            self.axes.set_title('Feret diameter of ' + type_figure)
            self.axes.legend()
            self.draw()

        if characteristic =='aspect ratio':
            self.axes.plot(*extract_and_trim_dict(aspect))
            self.axes.set_xlabel('aspect ratio')
            self.axes.set_ylabel('p')
            self.axes.set_title('Aspect ratio of ' + type_figure)
            self.draw()

        if characteristic =='sphericity':
            self.axes.plot(*extract_and_trim_dict(sphericity))
            self.axes.set_xlabel('sphericity')
            self.axes.set_ylabel('p')
            self.axes.set_title('Sphericity of ' + type_figure)
            self.draw()
                
                

                
class WidgetGallery(QDialog):
    def __init__(self, parent=None):
        super(WidgetGallery, self).__init__(parent)

        self.originalPalette = QApplication.palette()


        self.createTopLeftGroupBox()
        self.createTopRightGroupBox()

        self.createBottomBox()
        self.createBottomRightBox()
        
        self.customCanvas = CustomCanvas(self)
        self.customCanvas.setGeometry(0, 0, 500, 500)
        self.customCanvas.setFixedSize(500, 500)
        
        self.progress_bar = QProgressBar(self)
        self.progress_bar.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.worker_thread = WorkerThread()
        self.worker_thread.update_signal.connect(self.update_progress)
        

        mainLayout = QGridLayout()
        mainLayout.addWidget(self.topLeftGroupBox, 0, 0)
        mainLayout.addWidget(self.topRightGroupBox, 0, 1)
        mainLayout.addWidget(self.bottomBox, 1, 0)  # Add the bottom box to the layout
        mainLayout.addWidget(self.bottomRightBox, 1, 1)
        
        mainLayout.addWidget(self.customCanvas, 0, 2, 3, 1, alignment=Qt.AlignmentFlag.AlignCenter)
        mainLayout.addWidget(self.progress_bar, 2, 0, 1, 2)
        
        mainLayout.setRowStretch(0, 1)
        mainLayout.setRowStretch(1, 1)
        mainLayout.setRowStretch(2, 1)
        mainLayout.setColumnStretch(0, 1)
        mainLayout.setColumnStretch(1, 1)
        mainLayout.setColumnStretch(2, 1)
        self.setLayout(mainLayout)
        self.setWindowTitle("Calculate distribution densities")




    def createTopLeftGroupBox(self):
        self.topLeftGroupBox = QGroupBox("Polyhedron")

        self.radioButton_parallelepiped = QRadioButton("Parallelepiped")
        self.radioButton_triprism = QRadioButton("Triangular prism")
        self.radioButton_hex = QRadioButton("Hexagonal prism")
        self.radioButton_parallelepiped.setChecked(True)
        
        #если пользователь нажал - меняем параметры
        self.radioButton_parallelepiped.clicked.connect(self.updateParams)
        self.radioButton_triprism.clicked.connect(self.updateParams)
        self.radioButton_hex.clicked.connect(self.updateParams)


        layout = QVBoxLayout()
        layout.addWidget(self.radioButton_parallelepiped)
        layout.addWidget(self.radioButton_triprism)
        layout.addWidget(self.radioButton_hex)
        layout.addStretch(1)
        self.topLeftGroupBox.setLayout(layout)
        
        
        
    def createTopRightGroupBox(self):
        self.topRightGroupBox = QGroupBox("Characteristic")

        self.radioButton_angle = QRadioButton("angle")
        self.radioButton_area = QRadioButton("area")
        self.radioButton_perimeter = QRadioButton("perimeter")
        self.radioButton_semiaxes = QRadioButton("semiaxes")
        self.radioButton_feret = QRadioButton("Feret diameter")
        self.radioButton_aspect = QRadioButton("aspect ratio")
        self.radioButton_sphericity = QRadioButton("sphericity")
        self.radioButton_angle.setChecked(True)
        
        self.radioButton_angle.clicked.connect(self.updateParams)
        self.radioButton_area.clicked.connect(self.updateParams)
        self.radioButton_perimeter.clicked.connect(self.updateParams)
        self.radioButton_semiaxes.clicked.connect(self.updateParams)
        self.radioButton_feret.clicked.connect(self.updateParams)
        self.radioButton_aspect.clicked.connect(self.updateParams)
        self.radioButton_sphericity.clicked.connect(self.updateParams)

        layout = QVBoxLayout()
        layout.addWidget(self.radioButton_angle)
        layout.addWidget(self.radioButton_area)
        layout.addWidget(self.radioButton_perimeter)
        layout.addWidget(self.radioButton_semiaxes)
        layout.addWidget(self.radioButton_feret)
        layout.addWidget(self.radioButton_aspect)
        layout.addWidget(self.radioButton_sphericity)

        layout.addStretch(1)
        self.topRightGroupBox.setLayout(layout)
        
    
    def createBottomBox(self):
        self.bottomBox = QGroupBox("Parameters")
        self.bottomBoxLayout = QGridLayout()
        
        self.param_n = QComboBox()
        self.param_n.addItems(['100', '1000', '10000', '100000', '1000000'])
        self.param_n.setCurrentIndex(0)
        self.param_a = QDoubleSpinBox()
        self.param_b = QDoubleSpinBox()
        self.param_c = QDoubleSpinBox()
        self.param_h = QDoubleSpinBox()
        self.param_r = QDoubleSpinBox()
        self.param_k = QDoubleSpinBox()
        
        self.param_h.setVisible(False)
        self.param_r.setVisible(False)
        self.param_k.setVisible(False)
        
        self.param_n_name = QLabel("n=")
        self.param_a_name = QLabel("a=")
        self.param_b_name = QLabel("b=")
        self.param_c_name = QLabel("c=")
        self.param_h_name = QLabel("h=")
        self.param_r_name = QLabel("r=")
        self.param_k_name = QLabel("k=")
        
        self.param_h_name.setVisible(False)
        self.param_r_name.setVisible(False)
        self.param_k_name.setVisible(False)
        
        
        self.param_a.setValue(1.)
        self.param_a.setMinimum(0.01)
        self.param_a.setMaximum(10.0)
        self.param_b.setValue(2.)
        self.param_b.setMinimum(0.01)
        self.param_b.setMaximum(10.0)
        self.param_c.setValue(3.)
        self.param_c.setMinimum(0.01)
        self.param_c.setMaximum(10.0)

        self.bottomBoxLayout.addWidget(self.param_n_name, 0, 0)
        self.bottomBoxLayout.addWidget(self.param_n, 0, 1)
        self.bottomBoxLayout.addWidget(self.param_a_name, 1, 0)
        self.bottomBoxLayout.addWidget(self.param_a, 1, 1)
        self.bottomBoxLayout.addWidget(self.param_b_name, 2, 0)
        self.bottomBoxLayout.addWidget(self.param_b, 2, 1)
        self.bottomBoxLayout.addWidget(self.param_c_name, 3, 0)
        self.bottomBoxLayout.addWidget(self.param_c, 3, 1)
        self.bottomBoxLayout.addWidget(self.param_h_name, 4, 0)
        self.bottomBoxLayout.addWidget(self.param_h, 4, 1)
        self.bottomBoxLayout.addWidget(self.param_r_name, 5, 0)
        self.bottomBoxLayout.addWidget(self.param_r, 5, 1)
        self.bottomBoxLayout.addWidget(self.param_k_name, 6, 0)
        self.bottomBoxLayout.addWidget(self.param_k, 6, 1)

        self.bottomBox.setLayout(self.bottomBoxLayout)
        
        
    def createBottomRightBox(self):
        
        self.bottomRightBox = QGroupBox()
        self.bottomRightBoxLayout = QGridLayout()
        
        self.calculate_button = QPushButton("Calculate!")
        self.calculate_button.clicked.connect(self.calculateAndPlot)
        
        self.btn_change_directory = QPushButton("Choose the directory to save distributions")
        self.btn_change_directory.clicked.connect(self.change_directory)
        self.current_directory = QTextEdit()
        
        self.bottomRightBoxLayout.addWidget(self.btn_change_directory)
        self.bottomRightBoxLayout.addWidget(QLabel('Current directory:'))
        self.bottomRightBoxLayout.addWidget(self.current_directory)
        self.bottomRightBoxLayout.addWidget(self.calculate_button)
        
        self.bottomRightBox.setLayout(self.bottomRightBoxLayout)
        
        cur_dir = os.getcwd()

        try:
            os.mkdir(cur_dir + '\\densities')
        except:
            pass
        
        self.current_directory.setText(str(cur_dir + '\\densities'))
        
    

    def change_directory(self):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.FileMode.Directory)
        dlg.setOption(QFileDialog.Option.ShowDirsOnly)

        if dlg.exec():
            filenames = dlg.selectedFiles()
            if filenames:
                self.current_directory.setText(str(filenames[0]))
        
        
    def update_progress(self, value):
        self.progress_bar.setValue(value)    
    
    
    def updateBottomBox(self, type_figure):
        
        self.param_n.setVisible(True)        
        self.param_a.setVisible(True)
        self.param_b.setVisible(True)
        self.param_c.setVisible(True)
        self.param_h.setVisible(True)
        self.param_r.setVisible(True)
        self.param_k.setVisible(True)
        
        self.param_a.setValue(1.)
        self.param_a.setMinimum(0.01)
        self.param_a.setMaximum(10.0)
        self.param_b.setValue(2.)
        self.param_b.setMinimum(0.01)
        self.param_b.setMaximum(10.0)
        self.param_c.setValue(3.)
        self.param_c.setMinimum(0.01)
        self.param_c.setMaximum(10.0)
        
        self.param_h.setValue(1.)
        self.param_h.setMinimum(0.01)
        self.param_h.setMaximum(10.0)
        
        self.param_r.setValue(1.) # 0 corresponds to a triangle at the base, 1 corresponds to a hexagon
        self.param_r.setMinimum(0.0)
        self.param_r.setMaximum(1.)
        self.param_k.setValue(1.) # specifies the height of the prism
        self.param_k.setMinimum(0.01)
        self.param_k.setMaximum(10.0)

        self.param_n_name.setVisible(True)        
        self.param_a_name.setVisible(True)
        self.param_b_name.setVisible(True)
        self.param_c_name.setVisible(True)
        self.param_h_name.setVisible(True)
        self.param_r_name.setVisible(True)
        self.param_k_name.setVisible(True)
        
        if type_figure=='parallelepiped':
        
            self.param_h.setVisible(False)
            self.param_r.setVisible(False)
            self.param_k.setVisible(False)
            
            self.param_h_name.setVisible(False)
            self.param_r_name.setVisible(False)
            self.param_k_name.setVisible(False)

        if type_figure=='triangular prism':
        
            self.param_a.setVisible(False)
            self.param_b.setVisible(False)
            self.param_c.setVisible(False)
            self.param_r.setVisible(False)
            self.param_k.setVisible(False)
            
            self.param_a_name.setVisible(False)
            self.param_b_name.setVisible(False)
            self.param_c_name.setVisible(False)
            self.param_r_name.setVisible(False)
            self.param_k_name.setVisible(False)
            
        if type_figure=='hex prism':
        
            self.param_a.setVisible(False)
            self.param_b.setVisible(False)
            self.param_c.setVisible(False)
            self.param_h.setVisible(False)
   
            self.param_a_name.setVisible(False)
            self.param_b_name.setVisible(False)
            self.param_c_name.setVisible(False)
            self.param_h_name.setVisible(False)

    

    
    def updateParams(self):
        
        if self.radioButton_parallelepiped.isChecked():
            self.updateBottomBox(type_figure='parallelepiped')

        if self.radioButton_triprism.isChecked():
            self.updateBottomBox(type_figure='triangular prism')
            
        if self.radioButton_hex.isChecked():
            self.updateBottomBox(type_figure='hex prism')
            
            
    def calculateAndPlot(self):
        
        path = self.current_directory.toPlainText()
        n = int(self.param_n.currentText())
        
        #Parallelepiped
        if self.radioButton_parallelepiped.isChecked() and self.radioButton_angle.isChecked():
            
            params = {'side_a': self.param_a.value(),
                      'side_b': self.param_b.value(),
                      'side_c': self.param_c.value()}
            
            #self.worker_thread.start() #start progress bar
            
            self.customCanvas.plot_distribution(type_figure='parallelepiped', 
                                                characteristic='angle', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_parallelepiped.isChecked() and self.radioButton_area.isChecked():
            
            params = {'side_a': self.param_a.value(),
                      'side_b': self.param_b.value(),
                      'side_c': self.param_c.value()}
            
            self.customCanvas.plot_distribution(type_figure='parallelepiped', 
                                                characteristic='area', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_parallelepiped.isChecked() and self.radioButton_perimeter.isChecked():
            
            params = {'side_a': self.param_a.value(),
                      'side_b': self.param_b.value(),
                      'side_c': self.param_c.value()}
            
            self.customCanvas.plot_distribution(type_figure='parallelepiped', 
                                                characteristic='perimeter', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)

        if self.radioButton_parallelepiped.isChecked() and self.radioButton_semiaxes.isChecked():
            
            params = {'side_a': self.param_a.value(),
                      'side_b': self.param_b.value(),
                      'side_c': self.param_c.value()}
            
            #self.worker_thread.start() #start progress bar
            
            self.customCanvas.plot_distribution(type_figure='parallelepiped', 
                                                characteristic='semiaxes', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_parallelepiped.isChecked() and self.radioButton_feret.isChecked():
            
            params = {'side_a': self.param_a.value(),
                      'side_b': self.param_b.value(),
                      'side_c': self.param_c.value()}
            
            self.customCanvas.plot_distribution(type_figure='parallelepiped', 
                                                characteristic='Feret diameter', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_parallelepiped.isChecked() and self.radioButton_aspect.isChecked():
            
            params = {'side_a': self.param_a.value(),
                      'side_b': self.param_b.value(),
                      'side_c': self.param_c.value()}
            
            self.customCanvas.plot_distribution(type_figure='parallelepiped', 
                                                characteristic='aspect ratio', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params) 
        
        if self.radioButton_parallelepiped.isChecked() and self.radioButton_sphericity.isChecked():
            
            params = {'side_a': self.param_a.value(),
                      'side_b': self.param_b.value(),
                      'side_c': self.param_c.value()}
            
            self.customCanvas.plot_distribution(type_figure='parallelepiped', 
                                                characteristic='sphericity', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params) 


        
        
        #Triangular prism
        if self.radioButton_triprism.isChecked() and self.radioButton_angle.isChecked():
            
            params = {'side': self.param_h.value()}
            
            self.customCanvas.plot_distribution(type_figure='triangular prism', 
                                                characteristic='angle', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_triprism.isChecked() and self.radioButton_area.isChecked():
            
            params = {'side': self.param_h.value()}
            
            self.customCanvas.plot_distribution(type_figure='triangular prism', 
                                                characteristic='area', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_triprism.isChecked() and self.radioButton_perimeter.isChecked():
            
            params = {'side': self.param_h.value()}
            
            self.customCanvas.plot_distribution(type_figure='triangular prism', 
                                                characteristic='perimeter', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_triprism.isChecked() and self.radioButton_semiaxes.isChecked():
            
            params = {'side': self.param_h.value()}
            
            self.customCanvas.plot_distribution(type_figure='triangular prism', 
                                                characteristic='semiaxes', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_triprism.isChecked() and self.radioButton_feret.isChecked():
            
            params = {'side': self.param_h.value()}
            
            self.customCanvas.plot_distribution(type_figure='triangular prism', 
                                                characteristic='Feret diameter', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_triprism.isChecked() and self.radioButton_aspect.isChecked():
            
            params = {'side': self.param_h.value()}
            
            self.customCanvas.plot_distribution(type_figure='triangular prism', 
                                                characteristic='aspect ratio', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)

        if self.radioButton_triprism.isChecked() and self.radioButton_sphericity.isChecked():
            
            params = {'side': self.param_h.value()}
            
            self.customCanvas.plot_distribution(type_figure='triangular prism', 
                                                characteristic='sphericity', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
                        
        #Hexagonal prism
        if self.radioButton_hex.isChecked() and self.radioButton_angle.isChecked():
            
            params = {'r': self.param_r.value(),
                      'k': self.param_k.value()}
            
            self.customCanvas.plot_distribution(type_figure='hex prism', 
                                                characteristic='angle', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_hex.isChecked() and self.radioButton_area.isChecked():
            
            params = {'r': self.param_r.value(),
                      'k': self.param_k.value()}
            
            self.customCanvas.plot_distribution(type_figure='hex prism', 
                                                characteristic='area', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_hex.isChecked() and self.radioButton_perimeter.isChecked():
            
            params = {'r': self.param_r.value(),
                      'k': self.param_k.value()}
            
            self.customCanvas.plot_distribution(type_figure='hex prism', 
                                                characteristic='perimeter', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
        
        if self.radioButton_hex.isChecked() and self.radioButton_semiaxes.isChecked():
            
            params = {'r': self.param_r.value(),
                      'k': self.param_k.value()}
            
            self.customCanvas.plot_distribution(type_figure='hex prism', 
                                                characteristic='semiaxes', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_hex.isChecked() and self.radioButton_feret.isChecked():
            
            params = {'r': self.param_r.value(),
                      'k': self.param_k.value()}
            
            self.customCanvas.plot_distribution(type_figure='hex prism', 
                                                characteristic='Feret diameter', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)
            
        if self.radioButton_hex.isChecked() and self.radioButton_aspect.isChecked():
            
            params = {'r': self.param_r.value(),
                      'k': self.param_k.value()}
            
            self.customCanvas.plot_distribution(type_figure='hex prism', 
                                                characteristic='aspect ratio', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)

        if self.radioButton_hex.isChecked() and self.radioButton_sphericity.isChecked():
            
            params = {'r': self.param_r.value(),
                      'k': self.param_k.value()}
            
            self.customCanvas.plot_distribution(type_figure='hex prism', 
                                                characteristic='sphericity', 
                                                PATH=path,
                                                N=n,
                                                thread=self.worker_thread,
                                                params=params)


if __name__ == '__main__':

    import sys

    app = QApplication(sys.argv)
    gallery = WidgetGallery()
    gallery.show()
    sys.exit(app.exec())                