from random import random, uniform, seed
from matplotlib import pyplot as plt
import numpy as np
import math
import os
import pickle
from scipy.spatial import ConvexHull

from .geom import distance, intersection, parallel, angle, orthogonal, solve
from .geom import Point, Vector, Line, Plane


def area_of_triangle_3d(triangle):
    """
    Функция для вычисления площади треугольника в трехмерном пространстве по координатам его вершин.

    Аргументы:
    v1, v2, v3 -- Point or np.darray, содержащие координаты вершин треугольника в трехмерном пространстве.

    Возвращает:
    Площадь треугольника.
    """
    v1 = triangle[0]
    v2 = triangle[1]
    v3 = triangle[2]
    if (isinstance(v1, Line) or isinstance(v2, Line) or isinstance(v3, Line)):
        return 100000
    
    # Вычисляем векторы, соединяющие вершины треугольника.
    a = (v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2])
    b = (v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2])
    # Вычисляем векторное произведение векторов a и b.
    cross_product = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    # Вычисляем длину вектора cross_product.
    length = (cross_product[0] ** 2 + cross_product[1] ** 2 + cross_product[2] ** 2)**0.5
    # Вычисляем площадь треугольника.
    area = length / 2
    return area

def area_of_polygon_3d(polygon):
    '''
    get area of polygon
    
    polygon: np.array of coordinates of points of polygon
    
    return: area of polygon
    '''
    #make triangulation
    triangles = from_polygon_to_list_of_triangles(polygon)
    S = 0.
    #for every triangle add area of polygon
    for triangle in triangles:
        S += area_of_triangle_3d(triangle)
    return S


def from_polygon_to_list_of_triangles(vectors):
    '''
    make triangulation of polygon
    
    vectors: np.array of coordinates of points of polygon
    
    return np.array of triangles
    '''
    n = len(vectors)
    fp = vectors[0]
    triangles = []
    for i in range(2, n):
        triangle = np.array([fp, vectors[i-1], vectors[i]])
        triangles.append(triangle)
    return np.array(triangles)


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

right_triangle = [Point(0, 0, 0), Point(1, 0, 0), Point(0.5, 3**0.5/2, 0)]

def is_point_in_triangle(point, polygon=right_triangle):
    '''
    Проверка находится ли точка в треугольнике в 3D
    '''
    #проверка находится ли точка внутри треугольника
    #проверяем равенство площадей
    S = area_of_triangle_3d(polygon)
    pA = polygon[0]
    pB = polygon[1]
    pC = polygon[2]
    S1 = area_of_triangle_3d([pA, pB, point])
    S2 = area_of_triangle_3d([pA, pC, point])
    S3 = area_of_triangle_3d([pC, pB, point])
    if abs(S-S1-S2-S3) < 0.0001:
        return True
    else:
        return False
    
    


def area2(i, j, k, poly):
    x1, y1 = poly[i][0], poly[i][1]
    x2, y2 = poly[j][0], poly[j][1]
    x3, y3 = poly[k][0], poly[k][1]
    
    area = 0.5 * abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)))
    return area


def dist2(i, j, poly):
    x1, y1 = poly[i][0], poly[i][1]
    x2, y2 = poly[j][0], poly[j][1]
    
    return ((x1-x2)**2 + (y1-y2)**2)**0.5


def next_point(p, poly):
    return (p + 1) % len(poly)


def antipodal_pairs(poly):
    '''
    assumed that poly is 2D convex
    '''
    ret = []
    ret_ix = []
    p = len(poly) - 1
    q = 0

    while (area2(p, next_point(p, poly), next_point(q, poly), poly) 
           > area2(p, next_point(p, poly), q, poly)):
        q = next_point(q, poly)
    
    q0 = q
    p0 = 0

    while q != p0:
        if len(ret) > 1000:
        #strange case
            print('strange things 1')
            return np.array(ret)[:10], np.array(ret_ix)[:10]
        p = next_point(p, poly)
        ret.append([poly[p], poly[q]])
        ret_ix.append([p, q])
        while (area2(p, next_point(p, poly), next_point(q, poly), poly)
               > area2(p, next_point(p, poly), q, poly)):
            q = next_point(q, poly)
            if len(ret) > 1000:
            #strange case
                print('strange things 2')
                return np.array(ret)[:10], np.array(ret_ix)[:10]
            if (p != q0) or (q != p0):
                ret.append([poly[p], poly[q]])
                ret_ix.append([p, q])
            else:
                return np.array(ret), np.array(ret_ix)
        if math.isclose(area2(p, next_point(p, poly), next_point(q, poly), poly),
                         area2(p, next_point(p, poly), q, poly)):
            # parallel case
            if len(ret) > 1000:
            #strange case
                print('strange things 3')
                return np.array(ret)[:10], np.array(ret_ix)[:10]
            if (p != q0) or (q != len(poly) - 1):
                ret.append([poly[p], poly[next_point(q, poly)]])
                ret_ix.append([p, next_point(q, poly)])
            else:
                ret.append([poly[next_point(p, poly)], poly[q]])
                ret_ix.append([next_point(p, poly), q])
                return np.array(ret), np.array(ret_ix)
    return np.array(ret), np.array(ret_ix)    
    
    
def from_rk_to_abwt(r, k, m=1):
    b = m/(2*r+1)
    a = r*b
    w = 3**0.5/2*(m-a)
    t = k*w
    return a, b, w, t


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
    
def make_random_plane_into_prism(side=1):
    '''
    генерирует плоскость внутри усечённой призмы равномерно в пространстве
    '''
    #генерируем плоскость по вектору нормали и свигу вдоль этой нормали

    center = Point(0.5, np.sqrt(3)/6, side/2)
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))

    prism = RegularTriangularPrism(side=side)
    
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


def make_random_plane_into_tetra(side_a=1, side_b=2, side_c=3):
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

    prism = Parallelepiped(side_a=side_a, side_b=side_b, side_c=side_c)
    
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

    
def make_random_plane_into_wcco(r, k, m=1):
    '''
    генерирует плоскость внутри усечённой призмы равномерно в пространстве
    '''
    #генерируем плоскость по вектору нормали и свигу вдоль этой нормали
    _,_,_,t = from_rk_to_abwt(r=r, k=k, m=m)
    center = Point(0.5, np.sqrt(3)/6, t/2)
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))

    wcco = WCCoPrism(r=r, k=k, m=m)
    
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
        
        self.ABC = Plane(self.A, self.B, self.C)
        self.FEG = Plane(self.F, self.E, self.G)
        self.ABF = Plane(self.A, self.B, self.F)
        self.BCG = Plane(self.B, self.C, self.G)
        self.CDH = Plane(self.C, self.D, self.H)
        self.ADE = Plane(self.A, self.D, self.E)
        
        self.iterable_edges = [self.AB, self.BC, self.AD, self.DC, self.AE, self.BF,
                               self.CG, self.DH, self.EF, self.EH, self.FG, self.HG]
        
        self.iterable_planes = [self.ABC, self.FEG, self.ABF, self.BCG, self.CDH, self.ADE]
        self.iterable_facets = [[self.A, self.B, self.C, self.D],
                               [self.E, self.F, self.G, self.H],
                               [self.A, self.B, self.F, self.E],
                               [self.B, self.C, self.G, self.F],
                               [self.C, self.D, self.H, self.G],
                               [self.A, self.D, self.H, self.E]]


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
    def __init__(self, side=1.):
        

        self.side = side
        a = side
        
        self.A = Point(0., 0., 0.)
        self.B = Point(1., 0., 0.)
        self.C = Point(0.5, 3**0.5/2, 0.)
        self.D = Point(0., 0., a)
        self.E = Point(1., 0., a)
        self.F = Point(0.5, 3**0.5/2, a)
        
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
        
        self.iterable_planes = [self.ABD, self.ACD, self.CBE, self.ABC, self.DEF]
        self.iterable_facets = [[self.A, self.B, self.E, self.D],
                               [self.A, self.C, self.F, self.D],
                               [self.C, self.B, self.E, self.F],
                               [self.A, self.B, self.C],
                               [self.D, self.E, self.F]]       
        
        
    def intersection_with_plane(self, plane):
        boundaries = self.iterable_edges

        intersections = filter(None, map(lambda edge: intersection(edge, plane), boundaries))
        intersections = filter(lambda x: not isinstance(x, Line), intersections)
        intersections = list(set(intersections))

        # Filter out any out of bounds intersections
        def in_bounds(point):
            #how we check it:
            #we need to check the z-coordinate of point: it should be between 0 and a
            #AND
            #the point [point[0], point[1], 0] shoud be inside the right triangle [A, B, C]
            
            right_triangle = [self.A, self.B, self.C]
            return (
                is_point_in_triangle(point=[point[0], point[1], 0], polygon=right_triangle)
                and self.A.z <= point.z <= self.D.z
            )
        
        intersections = list(filter(in_bounds, intersections))
        #ok we got a list of intersections, but we also need to order it (in polygon array).
        #so until we have anything in intersection we will add point in polygon and pop it from intersection
        #we will add a point which have a common plane
        
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
                if uu_kostyl > 10000:
                    print('i had a problem with precision')
                    return []
            return polygon
        else:
            return []


class Parallelepiped():
    #параллелепипед с вершиной А в нуле координат и заданной длиной ребра, находится в первой полуплоскости для каждой пары осей
    def __init__(self, side_a=1, side_b=2, side_c=3):
        
        
        self.side_a = side_a
        self.side_b = side_b
        self.side_c = side_c
        
        self.A = Point(0, 0, 0)
        self.B = Point(side_a, 0, 0)
        self.C = Point(side_a, side_b, 0)
        self.D = Point(0, side_b, 0)
        self.E = Point(0, 0, side_c)
        self.F = Point(side_a, 0, side_c)
        self.G = Point(side_a, side_b, side_c)
        self.H = Point(0, side_b, side_c)
        
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
        
        self.ABC = Plane(self.A, self.B, self.C)
        self.FEG = Plane(self.F, self.E, self.G)
        self.ABF = Plane(self.A, self.B, self.F)
        self.BCG = Plane(self.B, self.C, self.G)
        self.CDH = Plane(self.C, self.D, self.H)
        self.ADE = Plane(self.A, self.D, self.E)
        
        self.iterable_planes = [self.ABC, self.FEG, self.ABF, self.BCG, self.CDH, self.ADE]
        self.iterable_facets = [[self.A, self.B, self.C, self.D],
                               [self.E, self.F, self.G, self.H],
                               [self.A, self.B, self.F, self.E],
                               [self.B, self.C, self.G, self.F],
                               [self.C, self.D, self.H, self.G],
                               [self.A, self.D, self.H, self.E]]
#        
        
    #better to check once more    
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
                    for plane in self.iterable_planes:
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
    def __init__(self, r, k, m=1):
        
        self.a, self.b, self.w, self.t = from_rk_to_abwt(r, k, m=m)
        
        
        self.A = Point(self.a, 0, 0)
        self.B = Point(self.a+self.b, 0, 0)
        self.C = Point(1.5*self.a+self.b, 3**0.5/2*self.a, 0)
        self.D = Point(0.5*m+self.a/2, self.w, 0)
        self.E = Point(0.5*m-self.a/2, self.w, 0)
        self.F = Point(0.5*self.a, 3**0.5/2*self.a, 0)
        
        self.G = Point(self.a, 0, self.t)
        self.H = Point(self.a+self.b, 0, self.t)
        self.I = Point(1.5*self.a+self.b, 3**0.5/2*self.a, self.t)
        self.J = Point(0.5*m+self.a/2, self.w, self.t)
        self.K = Point(0.5*m-self.a/2, self.w, self.t)
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
        
        self.iterable_planes = [self.ABC, self.GHI, self.ABH, self.BCH, 
                                self.CDJ, self.DEK, self.EFL, self.FAG]
        
        self.iterable_facets = [[self.A, self.B, self.C, self.D, self.E, self.F],
                               [self.G, self.H, self.I, self.J, self.K, self.L],
                               [self.A, self.B, self.H, self.G],
                               [self.B, self.C, self.I, self.H],
                               [self.C, self.D, self.J, self.I],
                               [self.D, self.E, self.K, self.J],
                               [self.E, self.F, self.L, self.K],
                               [self.A, self.F, self.L, self.G]] 
        
    def is_point_in_truncated_triangle(self, point):
        '''
        check if point in [A, B, C, D, E, F]
        '''
        polygon = [self.A, self.B, self.C, self.D, self.E, self.F]
        #проверяем равенство площадей
        S = area_of_polygon_3d(polygon)
        pA = polygon[0]
        pB = polygon[1]
        pC = polygon[2]
        pD = polygon[3]
        pE = polygon[4]
        pF = polygon[5]
        S1 = area_of_triangle_3d([pA, pB, point])
        S2 = area_of_triangle_3d([pB, pC, point])
        S3 = area_of_triangle_3d([pC, pD, point])
        S4 = area_of_triangle_3d([pD, pE, point])
        S5 = area_of_triangle_3d([pE, pF, point])
        S6 = area_of_triangle_3d([pF, pA, point])
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
            #how we check it:
            #we need to check the z-coordinate of point: it should be between 0 and a
            #AND
            #the point [point[0], point[1], 0] shoud be inside the right triangle [A, B, C]
            return (
                self.is_point_in_truncated_triangle(point=[point[0], point[1], 0])
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


def secant_polygon(polyhedron, way, params=None):
    '''
    генерирует плоскость внутри многогранника polyhedron способом way
    и возвращает список вершин секущего многогранника
    рекомендуется выбирать way, полностью замощающий polyhedron
    в polyhedron должен быть реализован метод intersection_with_plane
    '''
    if params:
        plane = way(**params)
    else:
        plane = way()
    # отбираем только те плоскости у которых непустое пересечение с многогранником
    points_intersection = polyhedron.intersection_with_plane(plane)
    return points_intersection


def plot_secant_polygon(type_figure='cube'):
    '''
    построение 3D графика фигуры type_figure и ее случайного сечения
    '''
    if type_figure=='cube':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        center_cube = (0.5, 0.5, 0.5)
        side_cube = 1.0
        half_side_cube = side_cube / 2
        vertices_cube = [
            (center_cube[0] - half_side_cube, center_cube[1] - half_side_cube, center_cube[2] - half_side_cube),
            (center_cube[0] - half_side_cube, center_cube[1] - half_side_cube, center_cube[2] + half_side_cube),
            (center_cube[0] - half_side_cube, center_cube[1] + half_side_cube, center_cube[2] - half_side_cube),
            (center_cube[0] - half_side_cube, center_cube[1] + half_side_cube, center_cube[2] + half_side_cube),
            (center_cube[0] + half_side_cube, center_cube[1] - half_side_cube, center_cube[2] - half_side_cube),
            (center_cube[0] + half_side_cube, center_cube[1] - half_side_cube, center_cube[2] + half_side_cube),
            (center_cube[0] + half_side_cube, center_cube[1] + half_side_cube, center_cube[2] - half_side_cube),
            (center_cube[0] + half_side_cube, center_cube[1] + half_side_cube, center_cube[2] + half_side_cube)
        ]

        # Соединения вершин для куба
        edges_cube = [
            (0, 1), (1, 3), (3, 2), (2, 0),
            (4, 5), (5, 7), (7, 6), (6, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)
        ]

        # Рисуем грани
        for edge in edges_cube:
            x = [vertices_cube[edge[0]][0], vertices_cube[edge[1]][0]]
            y = [vertices_cube[edge[0]][1], vertices_cube[edge[1]][1]]
            z = [vertices_cube[edge[0]][2], vertices_cube[edge[1]][2]]
            ax.plot(x, y, z, color='purple')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.title('Куб со стороной '+ str(side_cube) + ' и центром ' + str(center_cube))
        
        vertices = secant_polygon(Cube(), make_random_plane_into_cube)
        print('Area is: ', area_of_polygon_3d(vertices))

        # Добавление линий между вершинами
        for i in range(len(vertices)):
            x1, y1, z1 = vertices[i]
            x2, y2, z2 = vertices[(i + 1) % len(vertices)]
            ax.plot([x1, x2], [y1, y2], [z1, z2], color='green')

        # Добавление вершин
        for vertex in vertices:
            x, y, z = vertex
            ax.scatter(x, y, z, color='blue')
        plt.show()
        
        
def draw_ellipse(center, radii, rotation_matrix, num_points=100):
    '''
    Функция отрисовки эллипса
    '''
    # Parameters of the ellipse
    theta = np.linspace(0, 2 * np.pi, num_points)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

    # Create the ellipse without rotation
    x = radii[0] * cos_theta
    y = radii[1] * sin_theta

    # Apply rotation using the rotation matrix
    points = np.vstack((x, y))
    rotated_points = np.dot(rotation_matrix, points)

    # Transpose the rotated points to get x and y
    x_rotated, y_rotated = rotated_points

    # Offset by the center
    x_rotated += center[0]
    y_rotated += center[1]

    # Plot the ellipse
    plt.figure(figsize=(6, 6))
    plt.plot(x_rotated, y_rotated)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('2D Ellipse')
    plt.grid(True)
    plt.axis('equal')
    print('Semiaxes: {:.2f}, {:.2f}'.format(radii[0], radii[1]))
        

def norm_dict(char_dict, step):
    '''
    Нормирует характеристику, поданную словарём, чтобы получить из гистограммы плотность вероятности
    Например, на вход подается гистограмма площадей:
    (Nebo от 24.10)
    
    {0.00: 5,
     0.01: 10,
     ...}
    step = 0.01
    
    Нужно каждое значение из словаря поделить на норм. константу
    norm_const = step*(values.sum())
    
    '''
    norm_const = step*sum(char_dict.values())
    return {key: value/norm_const for key, value in char_dict.items()} #генератор списков


#code from https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py
def getMinVolEllipse(P=None, tolerance=0.0001):
    """ Find the minimum volume ellipsoid which holds all the points

    Based on work by Nima Moshtagh
    http://www.mathworks.com/matlabcentral/fileexchange/9542
    and also by looking at:
    http://cctbx.sourceforge.net/current/python/scitbx.math.minimum_covering_ellipsoid.html
    Which is based on the first reference anyway!

    Here, P is a numpy array of N dimensional points like this:
    P = [[x,y,z,...], <-- one point per line
         [x,y,z,...],
         [x,y,z,...]]

    Returns:
    (center, radii, rotation)

    """
    (N, d) = np.shape(P)
    d = float(d)

    # Q will be our working array
    Q = np.vstack([np.copy(P.T), np.ones(N)]) 
    QT = Q.T

    # initializations
    err = 1.0 + tolerance
    u = (1.0 / N) * np.ones(N)

    # Khachiyan Algorithm
    while err > tolerance:
        V = np.dot(Q, np.dot(np.diag(u), QT))
        M = np.diag(np.dot(QT , np.dot(np.linalg.inv(V), Q)))    # M the diagonal vector of an NxN matrix
        j = np.argmax(M)
        maximum = M[j]
        step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0))
        new_u = (1.0 - step_size) * u
        new_u[j] += step_size
        err = np.linalg.norm(new_u - u)
        u = new_u

    # center of the ellipse 
    center = np.dot(P.T, u)

    # the A matrix for the ellipse
    A = np.linalg.inv(
                   np.dot(P.T, np.dot(np.diag(u), P)) - 
                   np.array([[a * b for b in center] for a in center])
                   ) / d

    # Get the values we'd like to return
    U, s, rotation = np.linalg.svd(A)
    radii = 1.0/np.sqrt(s)

    return (center, radii, rotation)


def find_rotation_matrix(Ax, Ay, Bx, By):
    '''
    get rotation matrix from one plane, containing vectors Ax, Ay 
    to another plane, containing vectors Bx, By
    https://math.stackexchange.com/questions/1876615/rotation-matrix-from-plane-a-to-b
    
    it works!! science!! is a power!!
    '''
    #make Ax атв Ay identity and orthogonal
    Ax = Ax/np.linalg.norm(Ax)
    Az = (np.cross(Ax, Ay))/np.linalg.norm(np.cross(Ax, Ay))
    Ay = (np.cross(Ax, Az))/np.linalg.norm(np.cross(Ax, Az))
    
    Bx = Bx/np.linalg.norm(Bx)
    Bz = (np.cross(Bx, By))/np.linalg.norm(np.cross(Bx, By))
    By = (np.cross(Bx, Bz))/np.linalg.norm(np.cross(Bx, Bz))
    
    #print(np.dot(Ax, Ay), np.dot(Bx, By)) 0.0 0.0
    
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    
    A[:, 0] = Ax
    A[:, 1] = Ay
    A[:, 2] = Az
    
    B[:, 0] = Bx
    B[:, 1] = By
    B[:, 2] = Bz
    
    return np.dot(B, np.linalg.inv(A))


def rotate_points_to_Oxy(points):
    '''
    polygon : (N, 3)
    
    return rotated_polygon(N, 2), rotation_matix, z-coordinate
    '''
    oxy_x = np.array([1., 0., 0])
    oxy_y = np.array([0., 1., 0])
    
    polygon = np.array([np.array(list(i)) for i in points])

    rotation_matrix = find_rotation_matrix(polygon[1] - polygon[0], polygon[2] - polygon[0], oxy_x, oxy_y)
    
    rotated_polygon = np.array([np.dot(rotation_matrix, point) for point in polygon])
    
    z = rotated_polygon[0][2]
    
    return rotated_polygon[:, :2], rotation_matrix, z
    
    
def rotate_points_back_to_3D(rotated_polygon, rotation_matrix, z):
    '''
    polygon : (N, 2)
    
    return polygon : (N, 3)
    '''
    polygon = np.zeros((len(rotated_polygon), 3))
    polygon[:, :2] = rotated_polygon;
    polygon[:, 2] = z
    polygon = np.array([np.dot(rotation_matrix.T, point) for point in polygon])
    return polygon


def is_convex_polygon_2D(points):
    """
    Проверяет, является ли многоугольник выпуклым.
    
    Аргументы:
    points: список точек, представленных в виде кортежей (x, y), задающих вершины многоугольника.
    
    Возвращает:
    True, если многоугольник выпуклый, False в противном случае.
    """
    n = len(points)
    if n < 3:
        return False  # Многоугольник должен содержать минимум 3 вершины

    # Проверяем направление оборота для каждой тройки последовательных вершин
    # Если знак угла между последовательными отрезками меняется, многоугольник не выпуклый
    def cross_product(p1, p2, p3):
        return (p2[0] - p1[0]) * (p3[1] - p2[1]) - (p2[1] - p1[1]) * (p3[0] - p2[0])

    direction = None
    for i in range(n):
        p1 = points[i]
        p2 = points[(i + 1) % n]
        p3 = points[(i + 2) % n]
        cp = cross_product(p1, p2, p3)
        if cp == 0:
            continue  # Точки коллинеарны, продолжаем проверку
        if direction is None:
            direction = cp > 0
        elif direction != (cp > 0):
            return False  # Знак угла изменился, многоугольник не выпуклый
    return True


def semiaxes_from_polygon_3D(polygon):
    '''
    polygon should contain at leats 3 points
    '''
    polygonOxy, _, _ = rotate_points_to_Oxy(polygon)
    _, radii, _ = getMinVolEllipse(polygonOxy)
    return min(radii), max(radii)


def perimeter_3D(points):
    polygon = np.array([np.array(list(i)) for i in points])
    per = 0.
    n = len(points)
    for i in range(n):
        per += np.linalg.norm(polygon[(i+1) % n]-polygon[i])
    return per


def minmax_feret(polygon):
    _, ret_ix = antipodal_pairs(polygon)
    #max feret diameter
    max_feret = 0
    for i in range(len(ret_ix)):
        temp = dist2(ret_ix[i][0], ret_ix[i][1], polygon)
        if temp > max_feret:
            max_feret = temp
    
    i = 0
    min_feret = 1000000.
    ret_ix = np.vstack([ret_ix, ret_ix[:, ::-1]])

    while i < len(ret_ix)-1: 
        v1 = ret_ix[i][0]  #0
        v2 = ret_ix[i+1][0] #1

        while v1 == v2:
            temp = 2*area2(v1, ret_ix[i][1], ret_ix[i+1][1], polygon) \
                        / dist2(ret_ix[i][1], ret_ix[i+1][1], polygon)
            if temp < min_feret:
                min_feret = temp

            #count area of (v1, ret_ix[i][1], ret_ix[i+1][1]) divide by length from ret_ix[i][1] to ret_ix[i+1][1]
            #if it is less than min - > remember it 

            i += 1
            v1 = ret_ix[i][0]  #1
            v2 = ret_ix[(i+1) % len(ret_ix)][0] #1

        #not equal first -> go to the next v1
        i += 1
    
    return min_feret, max_feret


def characteristic_of_section(secant_polygon):
    '''
    возвращает список углов, площадь, полуоси
    '''
    #углы, площадь, полуоси
    if len(secant_polygon) > 2:
        angles = get_angles(secant_polygon)
        area = area_of_polygon_3d(secant_polygon)
        perim = perimeter_3D(secant_polygon)
        a_semi, b_semi = semiaxes_from_polygon_3D(secant_polygon)
        
    else:
        angles = []
        area = 0.00
        a_semi, b_semi = 0.0, 0.0
        perim = 0.
  
    return angles, area, perim, a_semi, b_semi



def calculate_all_distributions(polyhedron, way, type_figure='parallelepiped', n=10000, params=None):
    
    #PATH = os.path.join(os.getcwd(), 'densities')

    os.makedirs('densities', exist_ok=True)

    #углы

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
        #path_ = os.path.join(PATH, f"parallelepiped_a={a_side}b={b_side}c={c_side}")
        path_ = os.path.join("densities", f"parallelepiped_a={a_side}b={b_side}c={c_side}")

        try:
            # angles = np.loadtxt(path_ + r'\anlges_N=' + str(n) + '.txt')
            # areas = np.load(path_ + r'\areas_N=' + str(n) + '.npy', allow_pickle=True).item()
            # a_semi = np.load(path_ + r'\asemi_N=' + str(n) + '.npy', allow_pickle=True).item()
            # b_semi = np.load(path_ + r'\bsemi_N=' + str(n) + '.npy', allow_pickle=True).item()
            # perimeter = np.load(path_ + r'\perimeter_N=' + str(n) + '.npy', allow_pickle=True).item()
            # #semiaxes = np.load(path_ + r'\semiaxes_N=' + str(n) + '.npy')
            # #ADD LOAD 4 PARAMS
            # minferet = np.load(path_ + r'\minferet_N=' + str(n) + '.npy', allow_pickle=True).item()
            # maxferet = np.load(path_ + r'\maxferet_N=' + str(n) + '.npy', allow_pickle=True).item()
            # aspect = np.load(path_ + r'\aspect_N=' + str(n) + '.npy', allow_pickle=True).item()
            # sphericity = np.load(path_ + r'\sphericity_N=' + str(n) + '.npy', allow_pickle=True).item()
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

            os.makedirs(path_, exist_ok=True)

            # np.savetxt(path_ + r'\anlges_N=' + str(n) + '.txt', angles)
            # np.save(path_ + r'\areas_N=' + str(n) + '.npy', areas)
            # np.save(path_ + r'\perimeter_N=' + str(n) + '.npy', perimeter)
            # np.save(path_ + r'\asemi_N=' + str(n) + '.npy', a_semi)
            # np.save(path_ + r'\bsemi_N=' + str(n) + '.npy', b_semi)
            # #np.save(path_ + r'\\semiaxes_N=' + str(n) + '.npy', semiaxes)

            # np.save(path_ + r'\minferet_N=' + str(n) + '.npy', minferet)
            # np.save(path_ + r'\maxferet_N=' + str(n) + '.npy', maxferet)
            # np.save(path_ + r'\aspect_N=' + str(n) + '.npy', aspect)
            # np.save(path_ + r'\sphericity_N=' + str(n) + '.npy', sphericity)

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

        #path_ = PATH + r'\tri-prism_side=' + str(side)
        path_ = os.path.join("densities", f"tri-prism_side={side}")

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

            os.makedirs(path_, exist_ok=True)

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

        path_ = os.path.join("densities", f"hex-prism_r={r}k={k}")

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

            os.makedirs(path_, exist_ok=True)

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


def make_distributions(polyhedron, way, PATH, type_figure='parallelepiped', n=10000, params=None):
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

                os.makedirs(path_, exist_ok=True)

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

                os.makedirs(path_, exist_ok=True)

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

                os.makedirs(path_, exist_ok=True)

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




def generate_point_in_triangle(A1, A2, A3):
    '''
    A1, A2, A3: Point or array - 2D/3D coordinates of triangle, where we generate a Point uniformly
    '''
    z = np.random.uniform(0., 1., 2)
    l1 = z[0]**0.5
    l2 = z[1]
    
    A1_ = np.array([A1[0], A1[1], A1[2]])
    A2_ = np.array([A2[0], A2[1], A2[2]])
    A3_ = np.array([A3[0], A3[1], A3[2]])
    
    P_new = (1-l1)*A1_ + l1*(1-l2)*A2_ + l1*l2*A3_
    return Point(P_new)


def is_point_in_polygon(P, polygon):
    '''
    check if point P is in polygon
    
    P: Point or np.darray
    polygon: np.darray of coordinates of polygon vertices
    
    return True if point is in polygon, False otherwise
    '''
    S = area_of_polygon_3d(polygon) #area of polygon
    S_ = 0 #area of triangles PAB, where A, B - neighbour points of polygon
    
    N_edges = len(polygon)
    for i in range(N_edges-1):
        S_ += area_of_triangle_3d([P, polygon[i], polygon[i+1]])
        
    S_ += area_of_triangle_3d([P, polygon[N_edges-1], polygon[0]])
    
    return abs(S-S_) < 0.0000001


def from_list_points_to_array(points):
    '''
    transform list of Points to array of 3D coordinates
    '''
    vecs = np.zeros((len(points), 3))
    for i in range(len(points)):
        vecs[i] = points[i][0], points[i][1], points[i][2]
    return vecs
    

def proportion_S(triangles):
    '''
    get array of proportion of triangles over the triangle
    
    triangles: list of triangles
    '''
    n_triangles = triangles.shape[0]
    proportions = np.zeros(n_triangles)
    
    for j in range(n_triangles):
        proportions[j] = area_of_triangle_3d(triangles[j])

    return proportions/proportions.sum()


def generations_point_in_polygon(vectors):
    '''
    generate point in polygon
    
    '''
    
    triangles = from_polygon_to_list_of_triangles(vectors)
    proportions = proportion_S(triangles)
    proportions_cumsum = np.cumsum(proportions) #to get prob of put point in one triangle
    
    num = np.random.uniform(0, 1)
    
    k = 0 #number of triangle
    for i in range(len(proportions_cumsum)):
        if num <= proportions_cumsum[i]:
            k = i
            break
    
    point = generate_point_in_triangle(triangles[k][0], triangles[k][1], triangles[k][2])
    return point


def test_generate_linear_intersept(cube=Cube()):
    '''
    возвращает список точек-проекций на плоскость и ее выпуклую оболочку
    '''
    #генерируем вектор нормали
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))
    
    #берем плоскость с таким вектором нормали и проходящую через 0.
    plane = Plane(Point(0., 0., 0.), n)
    #в этой плоскости надо найти многоугольник-проекцию всех точек куба
    polygon_points_not_ordered = []
    for M in cube.iterable_points:
        line_MP = Line(M, n)
        #P - projection of M on plane
        P = intersection(line_MP, plane)
        polygon_points_not_ordered.append([P[0], P[1], P[2]])
    
    polygon_points_not_ordered = np.array(polygon_points_not_ordered)
    hull = ConvexHull(polygon_points_not_ordered[:, :2])
    return polygon_points_not_ordered, hull



def generate_linear_intersept(polyhedron=Cube()):
    '''
    generation of length of random linear intersept in polyhedron
    '''
    #генерируем вектор нормали
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))
    
    #берем плоскость с таким вектором нормали и проходящую через 0.
    plane = Plane(Point(0., 0., 0.), n)
    #в этой плоскости надо найти многоугольник-проекцию всех точек куба
    polygon_points_not_ordered = []
    for M in polyhedron.iterable_points: #проходимся по точкам куба
        line_MP = Line(M, n) #линия-нормаль - через точку куба, которую проецируем
        P = intersection(line_MP, plane) #P - projection of M on plane
        polygon_points_not_ordered.append([P[0], P[1], P[2]]) #добавим точку P в список точек-проекций
    
    polygon_points_not_ordered = np.array(polygon_points_not_ordered) 
    hull = ConvexHull(polygon_points_not_ordered[:, :2]) #берем выпуклую оболочку у проекции многоугольника на Oxy чтобы найти
                                                         #вершины, которые являются тем самым многоугольником-проекцией
    polygon = polygon_points_not_ordered[hull.vertices]  #оставляем только нужные вершины и в нужном порядке
    
    Pns = generations_point_in_polygon(polygon) #генерируем точку в этом многоугольнике-проекции
    linear_intersept_line = Line(Pns, n) # прямая, содержащая этот перехват, проходит через Pns и сонаправлена n
    
    two_points = [] #будем надеяться что будет две точки у концов отрезка
    #будем среди пересечений линии linear_intersept_line искать те, которые пересекают именно грани
    
    for i, pl in enumerate(polyhedron.iterable_planes):
        #check the intersection of generated line and a polyhedron's facets
        X = intersection(linear_intersept_line, pl)
        #need to check is the X inside the facet
        if is_point_in_polygon(X, from_list_points_to_array(polyhedron.iterable_facets[i])):
            two_points.append(X)
    if len(two_points) != 2:
        print('smth went wrong, not 2 points')
        return None
    P1 = two_points[0]
    P2 = two_points[1]
    lin = P1-P2
    
    return (P1, P2), (lin[0]**2 + lin[1]**2 + lin[2]**2)**0.5


def calculate_linear_intercept(polyhedron, type_figure='parallelepiped', n=10000, params=None):
    
    #PATH = os.path.join(os.getcwd(), 'densities')

    try:
        os.mkdir('densities')
    except:
        pass


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

        path_ = os.path.join("densities", f"parallelepiped_a={a_side}b={b_side}c={c_side}")

        try:
            intercept = np.load(os.path.join(path_, f'intercept_N={n}.npy'), allow_pickle=True).item()


        except:
            intercept_dict = {} #пустой словать с ключами от 0.00 до a^2+b^2+c^2 #можно меньше
            step = 0.01
            current_key = 0.00
            while current_key <= 10*(a_side**2 + b_side**2 + c_side**2):
                intercept_dict.setdefault(round(current_key, 2), 0)
                current_key += step

            num_errors = 0   


            for i in range(n):

                _, temp = generate_linear_intersept(polyhedron=polyhedron)
                if temp is not None:
                    intercept_dict[round(temp, 2)] += 1
                else:
                    num_errors += 1

            print('precision is ', 1-num_errors/n)

            intercept = norm_dict(intercept_dict, step=0.01)

            try:
                os.mkdir(path_)
            except:
                pass

            np.save(os.path.join(path_, f'intercept_N={n}.npy'), intercept)



    if type_figure == 'triangular prism':

        if params:
            side = params['side']
        else:
            side = 1.
            params = {'side': 1.}

        #path_ = PATH + r'\tri-prism_side=' + str(side)
        path_ = os.path.join("densities", f"tri-prism_side={side}")

        try:
            intercept = np.load(os.path.join(path_, f'intercept_N={n}.npy'), allow_pickle=True).item()

        except:
            intercept_dict = {} #пустой словать с ключами от 0.00 до a^2+b^2+c^2 #можно меньше
            step = 0.01
            current_key = 0.00
            while current_key <= 10*(side*side):
                intercept_dict.setdefault(round(current_key, 2), 0)
                current_key += step

            num_errors = 0   


            for i in range(n):
                _, temp = generate_linear_intersept(polyhedron=polyhedron)
                if temp is not None:
                    intercept_dict[round(temp, 2)] += 1
                else:
                    num_errors += 1

            print('precision is ', 1-num_errors/n)

            intercept = norm_dict(intercept_dict, step=0.01)

            try:
                os.mkdir(path_)
            except:
                pass

            np.save(os.path.join(path_, f'intercept_N={n}.npy'), intercept)





    if type_figure == 'hex prism':

        r = params['r']
        k = params['k']

        path_ = os.path.join("densities", f"hex-prism_r={r}k={k}")

        try:
            intercept = np.load(os.path.join(path_, f'intercept_N={n}.npy'), allow_pickle=True).item()

        except:
            intercept_dict = {} #пустой словать с ключами от 0.00 до a^2+b^2+c^2 #можно меньше
            step = 0.01
            current_key = 0.00
            while current_key <= 10.:
                intercept_dict.setdefault(round(current_key, 2), 0)
                current_key += step

            num_errors = 0   


            for i in range(n):
                _, temp = generate_linear_intersept(polyhedron=polyhedron)
                if temp is not None:
                    intercept_dict[round(temp, 2)] += 1
                else:
                    num_errors += 1

            print('precision is ', 1-num_errors/n)

            intercept = norm_dict(intercept_dict, step=0.01)

            try:
                os.mkdir(path_)
            except:
                pass

            np.save(os.path.join(path_, f'intercept_N={n}.npy'), intercept)


    return intercept