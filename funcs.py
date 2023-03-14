from random import random, uniform, seed
from matplotlib import pyplot as plt
import numpy as np
import math
import pandas as pd

from numba import jit, njit, prange
import pickle

from calc import distance, intersection, parallel, angle, orthogonal
from line import Line
from plane import Plane
from point import Point
from solver import solve
from vector import Vector

from IPython.display import clear_output


__all__ = (
    "Line",
    "Plane",
    "Point",
    "Vector",
    "angle",
    "distance",
    "intersection",
    "orthogonal",
    "parallel",
    "solve",
)


def get_angles(points):
    # вычисление углов полигона
    cnt = points
    num_points = len(cnt)
    angles = []
    for i, point in enumerate(cnt):
        point1 = cnt[i - 1]
        point2 = cnt[i]
        point3 = cnt[(i + 1) % num_points]
        angles.append(int(Vector(point2, point1).angle(Vector(point2, point3))* 180 / np.pi))

    return np.array(angles)

def S_from_triangle(triangle):
    # по точкам треугольника получить площадь
    x1 = triangle[0][0]
    y1 = triangle[0][1]
    x2 = triangle[1][0]
    y2 = triangle[1][1]    
    x3 = triangle[2][0]
    y3 = triangle[2][1]
    return 0.5*abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))

right_triangle = [Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 3**0.5/2, 0)]

def is_point_in_triangle(point, polygon=right_triangle):
    #проверка находится ли точка внутри треугольника
    #проверяем равенство площадей
    S = S_from_triangle(polygon)
    pA = polygon[0]
    pB = polygon[1]
    pC = polygon[2]
    S1 = S_from_triangle([pA, pB, point])
    S2 = S_from_triangle([pA, pC, point])
    S3 = S_from_triangle([pC, pB, point])
    if abs(S-S1-S2-S3) < 0.0001:
        return True
    else:
        return False
    
def from_rk_to_abwt(r, k):
    b = 1/(2*r+1)
    a = r*b
    w = 3**0.5/2*(1-a)
    t = k*w
    return a, b, w, t


#для прямой на плоскости

def points_to_line(x1, y1, x2, y2):
    if x1 != x2:
        k = (y2-y1)/(x2-x1)
        b = y1-k*x1
    return k, b

@njit
def intersect_line_map(k, b, n=10, n_line=100):
    eps = 0.000001
    map_2d_now = np.zeros((n, n))
    '''
    if k>0:
        mini = -b/k
        maxi = (n - b)/k
    elif k<0:
        mini = (n - b)/k
        maxi = -b/k
    '''
    mini=0
    maxi=n
    x_iter = np.linspace(max(0, mini), min(maxi, n) - eps, n_line)
    for x_n in x_iter:
        y_n = k*x_n+b
        if 0 <= y_n < n:
            i, j = math.floor(x_n), math.floor(y_n)
            map_2d_now[i][j] = 1
    return map_2d_now

@njit
def make_2d_map(n=10, n_line=100, num_lines=100, norm=True, seed=3407):
    '''
    Генерируем угол и сдвиг, правильный способ
    n*n - число квадратиков разбиения
    n_line - на сколько частей разбивается сгенерированная прямая, когда по ней итерируются
    num_lines - число генерируемых прямых
    '''
    #фикс сид
    random.seed(seed)
    
    map_2d = np.zeros((n, n))
    for i in range(num_lines):
        k = np.tan(uniform(0, 2*np.pi))
        if k >= 0:
            b = uniform(-k*n, n)
        else:
            b = uniform(0, n-k*n)
        map_2d += intersect_line_map(k, b, n, n_line)
    
    if norm:
        print("std is ", np.std(map_2d)/map_2d.sum())
        return map_2d/map_2d.sum()
    else:
        return map_2d




#для шара

def make_random_plane_into_sphere(r=0.866, center=Point(0.5, 0.5, 0.5)):
    '''
    генерирует плоскость внутри шара равномерно в пространстве
    r - радиус шара
    center - центр шара
    '''
    #генерируем плоскость по вектору нормали и свигу вдоль этой нормали
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(r*np.sin(teta)*np.cos(fi), 
               r*np.sin(teta)*np.sin(fi),
               r*np.cos(teta))
    alpha = uniform(-1, 1) #сдвиг
    point = Point(alpha*n[0]+center[0], alpha*n[1]+center[1], alpha*n[2]+center[2])
    return Plane(point, n)


def make_random_plane_into_cube(side=1, center=Point(0.5, 0.5, 0.5)):
    '''
    генерирует плоскость внутри куба равномерно в пространстве
    side - сторона куба
    center - центр куба
    '''
    #генерируем плоскость по вектору нормали и свигу вдоль этой нормали
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))

    cube = Cube(side=side, center=center)
    
    dist_max = 0
    for M in cube.iterable_points:
        #хочу найти проекцию всех точек куба на сгенерированную прямую
        #эта прямая
        a = Line(center, n) #через точку и вектор
        #проводим через точку M плоскость, перпендикулярную этой прямой
        pl = Plane(M, n) #через точку и вектор нормали
        L = pl.intersection(a) #искомая точка - пересечение плоскости и сгенерированной прямой
        dist = distance(center, L)
        if dist > dist_max:
            dist_max = dist
            
    #alpha - сдвиг, его будем подбирать в зависимости от сгенерированной прямой
    alpha = uniform(-dist_max, dist_max)
    point = Point(alpha*n[0]+center[0], alpha*n[1]+center[1], alpha*n[2]+center[2])
    return Plane(point, n)
    # проверка показала, что при такой генерации мы действительно обязательно пересекаем куб!!
    # и что если чуть-чуть пошевелим, то уже не обязательно пересекаем куб
    
def make_random_plane_into_prism():
    '''
    генерирует плоскость внутри усечённой призмы равномерно в пространстве
    '''
    #генерируем плоскость по вектору нормали и свигу вдоль этой нормали

    center = Point(0.5, np.sqrt(3)/6, 1/2)
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))

    prism = RegularTriangularPrism()
    
    arr_projs =[]
    for M in prism.iterable_points:
        #хочу найти проекцию всех вершин призмы на сгенерированную прямую
        #эта прямая
        a = Line(center, n) #через точку и вектор
        #проводим через точку M плоскость, перпендикулярную этой прямой
        pl = Plane(M, n) #через точку и вектор нормали
        L = pl.intersection(a) #искомая точка - пересечение плоскости и сгенерированной прямой
        arr_projs.append(L)

    dist_max = 0
    M_r = center
    M_l = center
    #среди всех попарных точек найдем те, между которыми расстояние максимальное
    # это нужно, потому что такая призма не симметричная
    for M1 in arr_projs:
        for M2 in arr_projs:
            dist = distance(M1, M2)
            if dist > dist_max:
                dist_max = dist
                M_l = M1
                M_r = M2
            
    d_min = distance(M_l, center)
    d_max = distance(M_r, center)
    #меняем направление вектора нормали, чтобы он смотрел о M_l до M_r
    n = Vector((M_r[0]-M_l[0])/distance(M_r, M_l), (M_r[1]-M_l[1])/distance(M_r, M_l), (M_r[2]-M_l[2])/distance(M_r, M_l))
    
    #alpha - сдвиг, его будем подбирать в зависимости от сгенерированной прямой нормали
    alpha = uniform(-d_min, d_max)
    point = Point(alpha*n[0]+center[0], alpha*n[1]+center[1], alpha*n[2]+center[2])
    return Plane(point, n)
    #проверено, все ок
    
def make_random_plane_into_wcco(r, k):
    '''
    генерирует плоскость внутри усечённой призмы равномерно в пространстве
    '''
    #генерируем плоскость по вектору нормали и свигу вдоль этой нормали
    _,_,_,t = from_rk_to_abwt(r=r, k=k)
    center = Point(0.5, np.sqrt(3)/6, t/2)
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))

    wcco = WCCoPrism(r=r, k=k)
    
    arr_projs =[]
    for M in wcco.iterable_points:
        #хочу найти проекцию всех вершин призмы на сгенерированную прямую
        #эта прямая
        a = Line(center, n) #через точку и вектор
        #проводим через точку M плоскость, перпендикулярную этой прямой
        pl = Plane(M, n) #через точку и вектор нормали
        L = pl.intersection(a) #искомая точка - пересечение плоскости и сгенерированной прямой
        arr_projs.append(L)

    dist_max = 0
    M_r = center
    M_l = center
    #среди всех попарных точек найдем те, между которыми расстояние максимальное
    # это нужно, потому что такая призма не симметричная
    for M1 in arr_projs:
        for M2 in arr_projs:
            dist = distance(M1, M2)
            if dist > dist_max:
                dist_max = dist
                M_l = M1
                M_r = M2
            
    d_min = distance(M_l, center)
    d_max = distance(M_r, center)
    #меняем направление вектора нормали, чтобы он смотрел о M_l до M_r
    n = Vector((M_r[0]-M_l[0])/distance(M_r, M_l), (M_r[1]-M_l[1])/distance(M_r, M_l), (M_r[2]-M_l[2])/distance(M_r, M_l))
    
    #alpha - сдвиг, его будем подбирать в зависимости от сгенерированной прямой нормали
    alpha = uniform(-d_min, d_max)
    point = Point(alpha*n[0]+center[0], alpha*n[1]+center[1], alpha*n[2]+center[2])
    return Plane(point, n)
    


def generate_random_plane_angles(polyhedron, way):
    '''
    генерирует плоскость внутри многогранника polyhedron способом way
    и возвращает список углов секущего многогранника
    рекомендуется выбирать way, полностью замощающий polyhedron
    в polyhedron должен быть реализован метод intersection_with_plane
    '''
    plane = way()
    # отбираем только те плоскости у которых непустое пересечение с многогранником
    points_intersection = polyhedron.intersection_with_plane(plane)
    if len(points_intersection) != 0:
        return get_angles(points_intersection)
    else:
        return []


def generate_n_planes(polyhedron, way, n=10000):
    '''
    генерирует n плоскостей способом way и возвращает плотность распределения
    углов секущих многогранников
    '''
    angles = [0]*181
    
    for i in range(n):
        angs = generate_random_plane_angles(polyhedron=polyhedron, way=way)
        if len(angs) != 0:
            for ang in angs:
                angles[ang] += 1
    
    return list(np.array(angles)/np.sum(angles))


def generate_random_plane_angles_wcco(polyhedron, way, r, k):
    '''
    генерирует плоскость внутри многогранника polyhedron способом way
    и возвращает список углов секущего многогранника
    рекомендуется выбирать way, полностью замощающий polyhedron
    в polyhedron должен быть реализован метод intersection_with_plane
    '''
    plane = way(r, k)
    # отбираем только те плоскости у которых непустое пересечение с многогранником
    points_intersection = polyhedron.intersection_with_plane(plane)
    if len(points_intersection) != 0:
        return get_angles(points_intersection)
    else:
        return []


def generate_n_planes_wcco(polyhedron, way, r, k, n=10000):
    '''
    генерирует n плоскостей способом way и возвращает плотность распределения
    углов секущих многогранников
    '''
    angles = [0]*181
    
    for i in range(n):
        angs = generate_random_plane_angles_wcco(polyhedron=polyhedron, way=way, r=r, k=k)
        if len(angs) != 0:
            for ang in angs:
                angles[ang] += 1
    
    return list(np.array(angles)/np.sum(angles))


def generate_random_plane_angles_type_wcco(polyhedron, way, r, k):
    '''
    генерирует плоскость внутри многогранника polyhedron способом way
    и возвращает список углов секущего многогранника и тип пересечения 
    0 - не пересекает базальные плоскости
    1 - одно пересечение
    2 - два пересечения
    рекомендуется выбирать way, полностью замощающий polyhedron
    в polyhedron должен быть реализован метод intersection_with_basal
    '''
    plane = way(r, k)
    # отбираем только те плоскости у которых непустое пересечение с многогранником
    points_intersection, tpe = polyhedron.intersection_with_basal(plane)
    if len(points_intersection) != 0:
        return get_angles(points_intersection), tpe
    else:
        return [], tpe
    
    
def generate_n_planes_type_wcco(polyhedron, way, r, k, n=10000):
    '''
    генерирует n плоскостей способом way и возвращает плотность распределения
    углов секущих многогранников
    '''
    angles = np.zeros((3, 181)) #= [0]*181
    
    for i in range(n):
        angs, tpe = generate_random_plane_angles_type_wcco(polyhedron=polyhedron, way=way, r=r, k=k)
        if len(angs) != 0:
            for ang in angs:
                angles[tpe][ang] += 1
    
    return list(np.array(angles)/np.sum(angles))


class Cube():
    #куб с вершиной А в нуле координат и заданной длиной ребра, находится в первой полуплоскости для каждой пары осей
    def __init__(self, side=1, center=Point(0.5, 0.5, 0.5)):
        
        
        self.side = side
        self.center = center
        a = side
        
        vec = Point(center - Point(side/2, side/2, side/2))
        
        self.A = Point(0, 0, 0) + vec
        self.B = Point(a, 0, 0) + vec
        self.C = Point(a, a, 0) + vec
        self.D = Point(0, a, 0) + vec
        self.E = Point(0, 0, a) + vec
        self.F = Point(a, 0, a) + vec
        self.G = Point(a, a, a) + vec
        self.H = Point(0, a, a) + vec
        
        self.iterable_points = [self.A, self.B, self.C, self.D, self.E, self.F, self.G, self.H]
        self.names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        
        self.AB = Line(self.A, self.B)
        self.BC = Line(self.B, self.C)
        self.AD = Line(self.A, self.D)
        self.DC = Line(self.D, self.C)
        
        self.AE = Line(self.A, self.E)
        self.BF = Line(self.B, self.F)
        self.CG = Line(self.C, self.G)
        self.DH = Line(self.D, self.H)
        
        self.EF = Line(self.E, self.F)
        self.EH = Line(self.E, self.H)
        self.FG = Line(self.F, self.G)
        self.HG = Line(self.H, self.G)
        
        self.iterable_edges = [self.AB, self.BC, self.AD, self.DC, self.AE, self.BF,
                               self.CG, self.DH, self.EF, self.EH, self.FG, self.HG]
        
        
        
    def intersection_with_plane(self, plane):
        boundaries = self.iterable_edges

        intersections = filter(None, map(lambda edge: intersection(edge, plane), boundaries))
        intersections = filter(lambda x: not isinstance(x, Line), intersections)
        intersections = list(set(intersections))

        # Filter out any out of bounds intersections
        def in_bounds(point):
            # intersect is actually (num, point)
            return (
                # <3 Python's comparison operator
                self.A.x <= point.x <= self.G.x and
                self.A.y <= point.y <= self.G.y and
                self.A.z <= point.z <= self.G.z
            )
        intersections = list(filter(in_bounds, intersections))
        
        
        if intersections:
            polygon = [intersections.pop()]
            while intersections:
                last = polygon[-1]
                distances = [distance(last, x) for x in intersections]
                # We're only interested in the index of the next point,
                # this min function returns the minimum (index, distance)
                # tuple...
                successor = min(enumerate(distances), key=lambda x: x[1])
                # ...but we only need the index :)
                successor = successor[0]
                polygon.append(intersections.pop(successor))

            return polygon
        else:
            return []
        

class RegularTriangularPrism():
    #Правильная треугольная призма, 
    def __init__(self):
        
        side=1
        
        center=Point(0.5, 0.5, 0.5)
        self.side = side
        self.center = center
        a = side
        
        vec = Point(center - Point(side/2, side/2, side/2))
        
        self.A = Point(0, 0, 0) + vec
        self.B = Point(a, 0, 0) + vec
        self.C = Point(0.5*a, 3**0.5/2*a, 0) + vec
        self.D = Point(0, 0, a) + vec
        self.E = Point(a, 0, a) + vec
        self.F = Point(0.5*a, 3**0.5/2*a, a) + vec
        
        self.iterable_points = [self.A, self.B, self.C, self.D, self.E, self.F]
        self.names = ['A', 'B', 'C', 'D', 'E', 'F']
        
        self.AB = Line(self.A, self.B)
        self.BC = Line(self.B, self.C)
        self.AC = Line(self.A, self.C)
        
        self.AD = Line(self.A, self.D)
        self.BE = Line(self.B, self.E)
        self.CF = Line(self.C, self.F)
        
        self.DE = Line(self.D, self.E)
        self.DF = Line(self.D, self.F)
        self.FE = Line(self.F, self.E)
        
        self.ABD = Plane(self.A, self.B, self.E)
        self.ACD = Plane(self.A, self.C, self.D)
        self.CBE = Plane(self.C, self.B, self.E)
        
        self.ABC = Plane(self.A, self.B, self.C)
        self.DEF = Plane(self.D, self.E, self.F)
        

        
        self.iterable_edges = [self.AB, self.BC, self.AC, self.AD, self.BE, self.CF,
                               self.DE, self.DF, self.FE]
        
        self.iterable_facets = [self.ABD, self.ACD, self.CBE, self.ABC, self.DEF]
        
        
        
    def intersection_with_plane(self, plane):
        boundaries = self.iterable_edges

        intersections = filter(None, map(lambda edge: intersection(edge, plane), boundaries))
        intersections = filter(lambda x: not isinstance(x, Line), intersections)
        intersections = list(set(intersections))

        # Filter out any out of bounds intersections
        def in_bounds(point):
            right_triangle = [self.A, self.B, self.C]
            # intersect is actually (num, point)
            return (
                is_point_in_triangle(point=point, polygon=right_triangle)
                and self.A.z <= point.z <= self.D.z
            )
        
        intersections = list(filter(in_bounds, intersections))
        #print('problem not here')
        
        if intersections:
            #print(intersections)
            point = intersections[0]
            polygon = []
            polygon.append(point)
            del intersections[0]
            
            uu_kostyl = 0 

            while intersections:
                for i, p in enumerate(intersections):
                    uu_kostyl += 1
                    flag = False
                    for plane in self.iterable_facets:
                        if (p in plane) and (point in plane):
                            polygon.append(p)
                            del intersections[i]
                            point = p
                            flag = True
                            break
                    if flag:
                        break
                if uu_kostyl > 10000:
                    print('i had a problem with precision')
                    return []
            return polygon
        else:
            return []

class WCCoPrism():
    def __init__(self, r, k):
        
        self.a, self.b, self.w, self.t = from_rk_to_abwt(r, k)
        
        
        self.A = Point(self.a, 0, 0)
        self.B = Point(self.a+self.b, 0, 0)
        self.C = Point(1.5*self.a+self.b, 3**0.5/2*self.a, 0)
        self.D = Point(0.5+self.a/2, self.w, 0)
        self.E = Point(0.5-self.a/2, self.w, 0)
        self.F = Point(0.5*self.a, 3**0.5/2*self.a, 0)
        
        self.G = Point(self.a, 0, self.t)
        self.H = Point(self.a+self.b, 0, self.t)
        self.I = Point(1.5*self.a+self.b, 3**0.5/2*self.a, self.t)
        self.J = Point(0.5+self.a/2, self.w, self.t)
        self.K = Point(0.5-self.a/2, self.w, self.t)
        self.L = Point(0.5*self.a, 3**0.5/2*self.a, self.t)
        
        
        self.iterable_points = [self.A, self.B, self.C, self.D, self.E, self.F,
                               self.G, self.H, self.I, self.J, self.K, self.L]
        self.names = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']
        
        self.AB = Line(self.A, self.B)
        self.BC = Line(self.B, self.C)
        self.CD = Line(self.C, self.D)
        self.DE = Line(self.D, self.E)
        self.EF = Line(self.E, self.F)
        self.FA = Line(self.F, self.A)
        
        self.AG = Line(self.A, self.G)
        self.BH = Line(self.B, self.H)
        self.CI = Line(self.C, self.I) 
        self.DJ = Line(self.D, self.J)
        self.EK = Line(self.E, self.K)
        self.FL = Line(self.F, self.L) 
        
        self.GH = Line(self.G, self.H)
        self.HI = Line(self.H, self.I)
        self.IJ = Line(self.I, self.J)
        self.JK = Line(self.J, self.K)
        self.KL = Line(self.K, self.L)
        self.LG = Line(self.L, self.G)
        
        self.ABC = Plane(self.A, self.B, self.C)
        self.GHI = Plane(self.G, self.H, self.I)
        
        self.ABH = Plane(self.A, self.B, self.H)
        self.BCH = Plane(self.B, self.C, self.H)
        self.CDJ = Plane(self.C, self.D, self.J)
        self.DEK = Plane(self.D, self.E, self.K)
        self.EFL = Plane(self.E, self.F, self.L)
        self.FAG = Plane(self.F, self.A, self.G)
        

        
        self.iterable_edges = [self.AB, self.BC, self.CD, self.DE, self.EF, self.FA,
                               self.AG, self.BH, self.CI, self.DJ, self.EK, self.FL, 
                               self.GH, self.HI, self.IJ, self.JK, self.KL, self.LG]
        
        self.iterable_facets = [self.ABC, self.GHI, self.ABH, self.BCH, 
                                self.CDJ, self.DEK, self.EFL, self.FAG]
        
    def is_point_in_truncated_triangle(self, point):
        polygon = [self.A, self.B, self.C, self.D, self.E, self.F, self.G]
        #проверяем равенство площадей
        S = 3**0.5/4-3*3**0.5/4*self.a**2
        pA = polygon[0]
        pB = polygon[1]
        pC = polygon[2]
        pD = polygon[3]
        pE = polygon[4]
        pF = polygon[5]
        S1 = S_from_triangle([pA, pB, point])
        S2 = S_from_triangle([pB, pC, point])
        S3 = S_from_triangle([pC, pD, point])
        S4 = S_from_triangle([pD, pE, point])
        S5 = S_from_triangle([pE, pF, point])
        S6 = S_from_triangle([pF, pA, point])
        if abs(S-S1-S2-S3-S4-S5-S6) < 0.000001:
            return True
        else:
            return False
        
        
    def intersection_with_plane(self, plane):
        boundaries = self.iterable_edges

        intersections = filter(None, map(lambda edge: intersection(edge, plane), boundaries))
        intersections = filter(lambda x: not isinstance(x, Line), intersections)
        intersections = list(set(intersections))

        # Filter out any out of bounds intersections
        def in_bounds(point):
            #right_triangle = [self.A, self.B, self.C]
            # intersect is actually (num, point)
            return (
                self.is_point_in_truncated_triangle(point=point)
                and 0 <= point.z <= self.t
            )
        
        intersections = list(filter(in_bounds, intersections))
        
        if intersections:
            point = intersections[0]
            polygon = []
            polygon.append(point)
            del intersections[0]
            
            uu_kostyl = 0 

            while intersections:
                for i, p in enumerate(intersections):
                    uu_kostyl += 1
                    flag = False
                    for plane in self.iterable_facets:
                        if (p in plane) and (point in plane):
                            polygon.append(p)
                            del intersections[i]
                            point = p
                            flag = True
                            break
                    if flag:
                        break
                if uu_kostyl > 100:
                    print('i had a problem with precision')
                    return []
            return polygon
        else:
            return []
        
    
    def intersection_with_basal(self, plane):
        polygon = self.intersection_with_plane(plane)
        if polygon:
            #basal planes
            #self.ABC, self.GHI
            polygon2 = polygon[1:] + [polygon[0]]
            flag_abc = False
            flag_ghi = False
            
            #пройдёмся по всем точкам многоугольника
            #и посмотрим, пересекают ли они базальные плоскости
            for i, point in enumerate(polygon):
                if Line(polygon[i], polygon2[i]) in self.ABC:
                    flag_abc = True
                if Line(polygon[i], polygon2[i]) in self.GHI:
                    flag_ghi = True
            
            if flag_abc and flag_ghi:
                return polygon, 2
            if flag_abc or flag_ghi:
                return polygon, 1
            else:
                return polygon, 0
        else:
            return polygon, -1


def generate_vertex_distribution(polyhedron, way, n=1000000):
    vertex = [0]*10
    for i in range(n):
        angs = generate_random_plane_angles(polyhedron, way)
        vertex[len(angs)] += 1
    return vertex

def generate_vertex_distribution_wcco(polyhedron, way, r, k, n=1000000):
    vertex = [0]*10
    for i in range(n):
        angs = generate_random_plane_angles_wcco(polyhedron, way, r, k)
        vertex[len(angs)] += 1
    return vertex

def make_table_vertex_distribution(polyhedron, way, n=100000, n_try=5):
    vertex_distrs = []
    for i in range(n_try):
        vertex_distrs.append(generate_vertex_distribution(polyhedron, way, n=n))
    vertex_distrs = np.array(vertex_distrs)
    vertex_distrs[:, 0] = 0
    vertex_distrs = vertex_distrs/np.sum(vertex_distrs, axis=1)[:,None]*100
    mean = np.mean(vertex_distrs, axis=0)
    std = np.std(vertex_distrs, axis=0)
    var = np.var(vertex_distrs, axis=0)
    
    indexes = list(range(1, n_try+1)) + ['mean', 'std', 'dispersion']
    return pd.DataFrame(np.vstack((vertex_distrs, mean, std, var)), index=indexes)

def make_table_vertex_distribution_wcco(polyhedron, way, r, k, n=100000, n_try=5):
    vertex_distrs = []
    for i in range(n_try):
        vertex_distrs.append(generate_vertex_distribution_wcco(polyhedron, way, r, k, n=n))
    vertex_distrs = np.array(vertex_distrs)
    vertex_distrs[:, 0] = 0
    vertex_distrs = vertex_distrs/np.sum(vertex_distrs, axis=1)[:,None]*100
    mean = np.mean(vertex_distrs, axis=0)
    std = np.std(vertex_distrs, axis=0)
    var = np.var(vertex_distrs, axis=0)
    
    indexes = list(range(1, n_try+1)) + ['mean', 'std', 'dispersion']
    return pd.DataFrame(np.vstack((vertex_distrs, mean, std, var)), index=indexes)


def get_path_from_r_k(r, k, n=100000):
    return "r=" + str(r) + "_k=" + str(k) + "_n=" + str(n)

def get_distribution_form_r_k(path, r, k, n=100000):
    '''
    посчитать распределение и записать на диск
    '''
    path_ = get_path_from_r_k(r, k, n)
    
    try:
        arr = np.loadtxt(path + "carbid_model_distributions/" + path_ + ".txt")
    except:
        arr = generate_n_planes_wcco(polyhedron=WCCoPrism(r=r, k=k), way=make_random_plane_into_wcco, r=r, k=k, n=n)
        
        path_to_save = path + "carbid_model_distributions/" + path_ + ".txt"
        np.savetxt(path_to_save, np.array(arr))
    return arr

def make_plot_from_r_k(path, r, k, n=100000):
    '''
    попытаться считать с диска, если для таких параметров уже посчитано распределение,
    то просто его выгрузить. Иначе посчитать.
    
    И сделать график
    '''
    
    path_ = get_path_from_r_k(r, k, n)
    description = "Распределение углов секущих многогранников для усечённой призмы с параметрами r=" + str(r) \
    + ", k=" + str(k)
    
    try:
        arr = np.loadtxt(path + "carbid_model_distributions/" + path_ + ".txt")
    except:
        arr = np.array(get_distribution_form_r_k(path, r=r, k=k, n=n))
    
    plt.figure(figsize=(7,7))
    plt.plot(arr)
    plt.title(description)
    plt.savefig(path + "carbid_model_plots/" + path_ + ".pdf", format='pdf')
    plt.show()
    #plt.clean_output()
    
    #0 1    99 [100 101 102    119 120 121 123] это индексы считала
    
    peaks = {'fa': arr[:100].argmax(), #first argument(about 90)
             'fm': arr[:100].max(),    #first max
             'sa': 100 + arr[100:].argmax(), #second argument(about 120)
             'sm': arr[100:].max()     #second max
            }
    
    #запишем в словать
    with open(path + 'peaks_dictionary/' + path_ + '.pkl', 'wb') as f:
        pickle.dump(peaks, f)

#     открыть словарь с пиками
#     with open('peaks_dictionary/' + path_ + '.pkl', 'rb') as f:
#         loaded_dict = pickle.load(f)
    
    return peaks


#загрузка всех пиков в один словарь для дальнейшей работы

def get_dict_all_peaks(path,
                       r_list = [0.2, 0.3, 0.5, 0.6, 0.7, 1],
                       k_list = [0.6, 0.75, 1, 1.5],
                       n = 100000):
    dict_all_peaks = {}
    for r in r_list:
        for k in k_list:
            path_ = get_path_from_r_k(r, k, n)
            with open(path + 'peaks_dictionary/' + path_ + '.pkl', 'rb') as f:
                dict_all_peaks[path_] = pickle.load(f)
    return dict_all_peaks