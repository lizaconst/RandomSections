from random import uniform
import numpy as np
import math
import os
from scipy.spatial import ConvexHull

from .geom import distance, intersection
from .geom import Point, Vector, Line, Plane


def area_of_triangle_3d(triangle):
    """
    Calculate the area of a triangle in 3D space given the coordinates of its vertices.

    Args:
        triangle (list or tuple of lists/tuples): A list or tuple containing three points, each of which is a list or tuple 
        representing the coordinates of a vertex of the triangle in 3D space.

    Returns:
        float: The area of the triangle.

    Notes:
        If any of the vertices is an instance of the `Line` class, the function returns a default value of 100000.

    """
    v1 = triangle[0]
    v2 = triangle[1]
    v3 = triangle[2]
    if (isinstance(v1, Line) or isinstance(v2, Line) or isinstance(v3, Line)):
        return 100000
    
    # Calculate vectors connecting the vertices of the triangle.
    a = (v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2])
    b = (v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2])
    # Compute the cross product of vectors a and b.
    cross_product = (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    # Calculate the length of the cross product vector.
    length = (cross_product[0] ** 2 + cross_product[1] ** 2 + cross_product[2] ** 2)**0.5
    area = length / 2
    return area


def area_of_polygon_3d(polygon):
    """
    Calculate the area of a polygon in 3D space given the coordinates of its vertices.

    Args:
        polygon (np.array): An array of coordinates representing the vertices of the polygon in 3D space.

    Returns:
        float: The area of the polygon.

    """
    # Make triangulation
    triangles = from_polygon_to_list_of_triangles(polygon)
    S = 0.0
    # For every triangle, add the area to the total polygon area
    for triangle in triangles:
        S += area_of_triangle_3d(triangle)
    return S


def from_polygon_to_list_of_triangles(vectors):
    """
    Make triangulation of a polygon.

    Args:
        vectors (np.array): An array of coordinates representing the vertices of the polygon in 3D space.

    Returns:
        np.array: An array of triangles, where each triangle is represented by the coordinates of its three vertices.
    """
    n = len(vectors)
    fp = vectors[0]
    triangles = []
    for i in range(2, n):
        triangle = np.array([fp, vectors[i-1], vectors[i]])
        triangles.append(triangle)
    return np.array(triangles)


def from_list_points_to_array(points):
    """
    Transform a list of Points to an array of 3D coordinates.

    Args:
        points (list): A list of Points, where each Point has 3D coordinates (x, y, z).

    Returns:
        np.array: An array of shape (n, 3) containing the 3D coordinates of the points.
    """
    vecs = np.zeros((len(points), 3))
    for i in range(len(points)):
        vecs[i] = points[i][0], points[i][1], points[i][2]
    return vecs


def from_list_array_to_points(points_arr):
    """
    Transform an array of 3D coordinates to a list of Points.

    Args:
        points_arr (np.array): An array of shape (n, 3) containing the 3D coordinates of the points.

    Returns:
        list: A list of Points, where each Point has 3D coordinates (x, y, z).
    """
    points = []
    for i in range(len(points_arr)):
        points.append(Point(points_arr[i][0], points_arr[i][1], points_arr[i][2]))
    return points


def get_angles(points):
    """
    Calculate the angles of a convex polygon at each vertex.

    Args:
        points (list or np.array or list (np.array) of points): A list or array of coordinates representing the vertices of the polygon in 3D space.

    Returns:
        np.array: An array of angles (in degrees) at each vertex of the polygon.

    """
    if isinstance(points[0], list) or isinstance(points[0], np.ndarray):
        points = from_list_array_to_points(points)

    cnt = points
    num_points = len(cnt)
    angles = []
    for i, point in enumerate(cnt):
        point1 = cnt[i - 1]
        point2 = cnt[i]
        point3 = cnt[(i + 1) % num_points]
        angles.append(int(Vector(point2, point1).angle(Vector(point2, point3))* 180 / np.pi))

    return np.array(angles)



def is_point_in_triangle(point, triangle):
    """
    Check if a point is inside a triangle in 3D space.

    Args:
        point (list or np.array): The coordinates of the point to be checked.
        triangle (list or np.array): An array of coordinates representing the vertices of the triangle in 3D space.

    Returns:
        bool: True if the point is inside the triangle, False otherwise.
    """

    S = area_of_triangle_3d(triangle)
    pA = triangle[0]
    pB = triangle[1]
    pC = triangle[2]
    S1 = area_of_triangle_3d([pA, pB, point])
    S2 = area_of_triangle_3d([pA, pC, point])
    S3 = area_of_triangle_3d([pC, pB, point])
    # Check if the sum of the sub-triangle areas is equal to the original triangle area
    if abs(S-S1-S2-S3) < 0.0001:
        return True
    else:
        return False
    

def antipodal_pairs(poly):
    """
    Find antipodal pairs in a 2D convex polygon.

    Args:
        poly (np.array): A 2D array representing the vertices of the convex polygon.

    Returns:
        tuple: Two numpy arrays:
            - An array of pairs of points representing the antipodal pairs.
            - An array of pairs of indices corresponding to the antipodal pairs.
    """

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
    """
    Convert parameters r and k of a truncated prism to parameters a, b, w, and t of a truncated prism.
    Parameter m is a scale parameter. By default, a + 2b = 1.

    Args:
        r (float): Parameter r of the truncated prism.
        k (float): Parameter k of the truncated prism.
        m (float, optional): Scale parameter. Defaults to 1.

    Returns:
        tuple: A tuple containing parameters a, b, w, and t of the truncated prism.

    """
    b = m/(2*r+1)
    a = r*b
    w = 3**0.5/2*(m-a)
    t = k*w
    return a, b, w, t


def extract_and_trim_dict(dictionary):
    """
    Trim a dictionary to remove trailing zeros and extract keys and values.

    Args:
        dictionary (dict): The dictionary to be trimmed.

    Returns:
        tuple: A tuple containing two lists, one with keys and the other with corresponding values after trimming.
    """
    keys = []
    values = []
    for key, value in dictionary.items():
        keys.append(key)
        values.append(value)

    while values and values[-1] == 0.0:
        del keys[-1]
        del values[-1]
    return keys, values


def make_random_plane_into_sphere(r=0.866, center=Point(0.5, 0.5, 0.5)):
    """
    Generate a plane inside a sphere uniformly in space.

    Args:
        r (float, optional): Radius of the sphere. Defaults to 0.866.
        center (Point, optional): Center of the sphere. Defaults to Point(0.5, 0.5, 0.5).

    Returns:
        Plane: A plane generated inside the sphere.

    """
    # Generate a plane along the normal vector and shift it along this normal vector
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(r*np.sin(teta)*np.cos(fi), 
               r*np.sin(teta)*np.sin(fi),
               r*np.cos(teta))
    alpha = uniform(-1, 1) #shift
    point = Point(alpha*n[0]+center[0], alpha*n[1]+center[1], alpha*n[2]+center[2])
    return Plane(point, n)


def make_random_plane_into_cube(side=1, center=Point(0.5, 0.5, 0.5)):
    """
    Generate a plane inside a cube uniformly in space.

    Args:
        side (float, optional): Side length of the cube. Defaults to 1.
        center (Point, optional): Center of the cube. Defaults to Point(0.5, 0.5, 0.5).

    Returns:
        Plane: A plane generated inside the cube.
    """
    # Generate a plane along the normal vector and shift it along this normal vector
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))

    cube = Cube(side=side, center=center)
    
    dist_max = 0
    for M in cube.iterable_points:
        # Find the projection of all points of the cube onto the generated line
        a = Line(center, n) # by point and vector
        #plane through M, perpendicular to the line n
        pl = Plane(M, n) #plane by point and normal vector
        L = pl.intersection(a) # intersection of plane and generated line
        dist = distance(center, L)
        if dist > dist_max:
            dist_max = dist
            
    # Choose the shift alpha depending on the generated line
    alpha = uniform(-dist_max, dist_max)
    point = Point(alpha*n[0]+center[0], alpha*n[1]+center[1], alpha*n[2]+center[2])
    return Plane(point, n)

    
def make_random_plane_into_prism(side=1):
    """
    Generate a plane inside a truncated prism uniformly in space.

    Args:
        side (float, optional): Side length of the truncated prism. Defaults to 1.

    Returns:
        Plane: A plane generated inside the truncated prism.
    """
    # Generate a plane along the normal vector and shift it along this normal vector
    center = Point(0.5, np.sqrt(3) / 6, side / 2)
    fi = 2 * np.pi * uniform(0, 1)
    teta = np.arccos(2 * uniform(0, 1) - 1)
    n = Vector(np.sin(teta) * np.cos(fi), 
               np.sin(teta) * np.sin(fi),
               np.cos(teta))

    prism = RegularTriangularPrism(side=side)
    
    arr_projs =[]
    for M in prism.iterable_points:
        # Find the projection of all vertices of the prism onto the generated line
        a = Line(center, n)
        pl = Plane(M, n)
        L = pl.intersection(a)
        arr_projs.append(L)

    dist_max = 0
    M_r = center
    M_l = center
    # Find the pair of points with maximum distance between them
    for M1 in arr_projs:
        for M2 in arr_projs:
            dist = distance(M1, M2)
            if dist > dist_max:
                dist_max = dist
                M_l = M1
                M_r = M2
            
    d_min = distance(M_l, center)
    d_max = distance(M_r, center)
    # Change the direction of the normal vector so that it points from M_l to M_r
    n = Vector((M_r[0] - M_l[0]) / distance(M_r, M_l), 
               (M_r[1] - M_l[1]) / distance(M_r, M_l), 
               (M_r[2] - M_l[2]) / distance(M_r, M_l))
    
    # Choose alpha depending on the generated normal line
    alpha = uniform(-d_min, d_max)
    point = Point(alpha * n[0] + center[0], 
                  alpha * n[1] + center[1], 
                  alpha * n[2] + center[2])
    return Plane(point, n)


def make_random_plane_into_tetra(side_a=1, side_b=2, side_c=3):
    """
    Generate a plane inside a parallelepiped uniformly in space.

    Args:
        side_a (float, optional): Length of side 'a' of the parallelepiped. Defaults to 1.
        side_b (float, optional): Length of side 'b' of the parallelepiped. Defaults to 2.
        side_c (float, optional): Length of side 'c' of the parallelepiped. Defaults to 3.

    Returns:
        Plane: A plane generated inside the parallelepiped.
    """
    # Generate a plane along the normal vector and shift it along this normal vector
    center = Point(0.5, np.sqrt(3) / 6, 1 / 2)
    fi = 2 * np.pi * uniform(0, 1)
    teta = np.arccos(2 * uniform(0, 1) - 1)
    n = Vector(np.sin(teta) * np.cos(fi), 
               np.sin(teta) * np.sin(fi),
               np.cos(teta))

    prism = Parallelepiped(side_a=side_a, side_b=side_b, side_c=side_c)
    
    arr_projs =[]
    for M in prism.iterable_points:
        # Find the projection of all vertices of the prism onto the generated line
        a = Line(center, n)
        pl = Plane(M, n)
        L = pl.intersection(a)
        arr_projs.append(L)

    dist_max = 0
    M_r = center
    M_l = center
    # Find the pair of points with maximum distance between them
    for M1 in arr_projs:
        for M2 in arr_projs:
            dist = distance(M1, M2)
            if dist > dist_max:
                dist_max = dist
                M_l = M1
                M_r = M2
            
    d_min = distance(M_l, center)
    d_max = distance(M_r, center)
    # Change the direction of the normal vector so that it points from M_l to M_r
    n = Vector((M_r[0] - M_l[0]) / distance(M_r, M_l), 
               (M_r[1] - M_l[1]) / distance(M_r, M_l), 
               (M_r[2] - M_l[2]) / distance(M_r, M_l))
    
    # Choose alpha depending on the generated normal line
    alpha = uniform(-d_min, d_max)
    point = Point(alpha * n[0] + center[0], 
                  alpha * n[1] + center[1], 
                  alpha * n[2] + center[2])
    return Plane(point, n)

    
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
    

class Cube():
    """
    Represents a cube with one vertex at the origin and a specified side length, located in the first octant for each pair of axes.

    Attributes:
        side (float): The length of the side of the cube.
        center (Point): The center point of the cube.
        A (Point): Vertex A of the cube.
        B (Point): Vertex B of the cube.
        C (Point): Vertex C of the cube.
        D (Point): Vertex D of the cube.
        E (Point): Vertex E of the cube.
        F (Point): Vertex F of the cube.
        G (Point): Vertex G of the cube.
        H (Point): Vertex H of the cube.
        iterable_points (list): List of all vertices of the cube.
        names (list): Names of the vertices of the cube.
        AB (Line): Line segment AB.
        BC (Line): Line segment BC.
        AD (Line): Line segment AD.
        DC (Line): Line segment DC.
        AE (Line): Line segment AE.
        BF (Line): Line segment BF.
        CG (Line): Line segment CG.
        DH (Line): Line segment DH.
        EF (Line): Line segment EF.
        EH (Line): Line segment EH.
        FG (Line): Line segment FG.
        HG (Line): Line segment HG.
        ABC (Plane): Plane ABC.
        FEG (Plane): Plane FEG.
        ABF (Plane): Plane ABF.
        BCG (Plane): Plane BCG.
        CDH (Plane): Plane CDH.
        ADE (Plane): Plane ADE.
        iterable_edges (list): List of all edges of the cube.
        iterable_planes (list): List of all planes of the cube.
        iterable_facets (list): List of all facets of the cube.

    Methods:
        intersection_with_plane(plane): Finds the intersection points of the cube with a given plane.
    """

    def __init__(self, side=1, center=Point(0.5, 0.5, 0.5)):

        """
        Initialize a Cube object with the specified side length and center point.

        Args:
            side (float, optional): The length of the side of the cube. Defaults to 1.
            center (Point, optional): The center point of the cube. Defaults to Point(0.5, 0.5, 0.5).
        """

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
                successor = min(enumerate(distances), key=lambda x: x[1])
                # but we only need the index
                successor = successor[0]
                polygon.append(intersections.pop(successor))

            return polygon
        else:
            return []
        

class RegularTriangularPrism():
    """
    Represents a regular triangular prism.

    Attributes:
        side (float): The length of the side of the triangular base.
        A (Point): Vertex A of the prism.
        B (Point): Vertex B of the prism.
        C (Point): Vertex C of the prism.
        D (Point): Vertex D of the prism.
        E (Point): Vertex E of the prism.
        F (Point): Vertex F of the prism.
        iterable_points (list): List of all vertices of the prism.
        names (list): Names of the vertices of the prism.
        AB (Line): Line segment AB.
        BC (Line): Line segment BC.
        AC (Line): Line segment AC.
        AD (Line): Line segment AD.
        BE (Line): Line segment BE.
        CF (Line): Line segment CF.
        DE (Line): Line segment DE.
        DF (Line): Line segment DF.
        FE (Line): Line segment FE.
        ABD (Plane): Plane ABD.
        ACD (Plane): Plane ACD.
        CBE (Plane): Plane CBE.
        ABC (Plane): Plane ABC.
        DEF (Plane): Plane DEF.
        iterable_edges (list): List of all edges of the prism.
        iterable_planes (list): List of all planes of the prism.
        iterable_facets (list): List of all facets of the prism.

    Methods:
        intersection_with_plane(plane): Finds the intersection points of the prism with a given plane.

    """
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
                is_point_in_triangle(point=[point[0], point[1], 0], triangle=right_triangle)
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
                    for plane in self.iterable_planes:
                        if (p in plane) and (point in plane):
                            polygon.append(p)
                            del intersections[i]
                            point = p
                            flag = True
                            break
                    if flag:
                        break
                if uu_kostyl > 1000:
                    #print('i had a problem with precision')
                    return []
            return polygon
        else:
            return []


class Parallelepiped():
    """
    Represents a parallelepiped.

    Attributes:
        side_a (float): The length of the parallelepiped along the x-axis.
        side_b (float): The length of the parallelepiped along the y-axis.
        side_c (float): The length of the parallelepiped along the z-axis.
        A (Point): Vertex A of the parallelepiped.
        B (Point): Vertex B of the parallelepiped.
        C (Point): Vertex C of the parallelepiped.
        D (Point): Vertex D of the parallelepiped.
        E (Point): Vertex E of the parallelepiped.
        F (Point): Vertex F of the parallelepiped.
        G (Point): Vertex G of the parallelepiped.
        H (Point): Vertex H of the parallelepiped.
        iterable_points (list): List of all vertices of the parallelepiped.
        names (list): Names of the vertices of the parallelepiped.
        AB (Line): Line segment AB.
        BC (Line): Line segment BC.
        AD (Line): Line segment AD.
        DC (Line): Line segment DC.
        AE (Line): Line segment AE.
        BF (Line): Line segment BF.
        CG (Line): Line segment CG.
        DH (Line): Line segment DH.
        EF (Line): Line segment EF.
        EH (Line): Line segment EH.
        FG (Line): Line segment FG.
        HG (Line): Line segment HG.
        iterable_edges (list): List of all edges of the parallelepiped.
        ABC (Plane): Plane ABC.
        FEG (Plane): Plane FEG.
        ABF (Plane): Plane ABF.
        BCG (Plane): Plane BCG.
        CDH (Plane): Plane CDH.
        ADE (Plane): Plane ADE.
        iterable_planes (list): List of all planes of the parallelepiped.
        iterable_facets (list): List of all facets of the parallelepiped.

    Methods:
        intersection_with_plane(plane): Finds the intersection points of the parallelepiped with a given plane.

    """
    def __init__(self, side_a=1, side_b=2, side_c=3):
        """
        Initialize a Parallelepiped object with the specified side lengths.

        Args:
            side_a (float, optional): The length of the parallelepiped along the x-axis. Defaults to 1.0.
            side_b (float, optional): The length of the parallelepiped along the y-axis. Defaults to 2.0.
            side_c (float, optional): The length of the parallelepiped along the z-axis. Defaults to 3.0.
        """
        
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
                if uu_kostyl > 1000:
                    #print('i had a problem with precision')
                    return []
            return polygon
        else:
            return []


class WCCoPrism():
    """
    Represents a truncated triangular prism.

    Attributes:
        a (float): Length of the first segment.
        b (float): Length of the second segment.
        w (float): Width of the truncated triangular prism.
        t (float): Height of the truncated triangular prism.
        A (Point): Vertex A of the prism base.
        B (Point): Vertex B of the prism base.
        C (Point): Vertex C of the prism base.
        D (Point): Vertex D of the prism base.
        E (Point): Vertex E of the prism base.
        F (Point): Vertex F of the prism base.
        G (Point): Vertex G of the top face.
        H (Point): Vertex H of the top face.
        I (Point): Vertex I of the top face.
        J (Point): Vertex J of the top face.
        K (Point): Vertex K of the top face.
        L (Point): Vertex L of the top face.
        iterable_points (list): List of all vertices of the prism.
        names (list): Names of the vertices of the prism.
        AB (Line): Line segment AB.
        BC (Line): Line segment BC.
        CD (Line): Line segment CD.
        DE (Line): Line segment DE.
        EF (Line): Line segment EF.
        FA (Line): Line segment FA.
        AG (Line): Line segment AG.
        BH (Line): Line segment BH.
        CI (Line): Line segment CI.
        DJ (Line): Line segment DJ.
        EK (Line): Line segment EK.
        FL (Line): Line segment FL.
        GH (Line): Line segment GH.
        HI (Line): Line segment HI.
        IJ (Line): Line segment IJ.
        JK (Line): Line segment JK.
        KL (Line): Line segment KL.
        LG (Line): Line segment LG.
        ABC (Plane): Plane ABC.
        GHI (Plane): Plane GHI.
        ABH (Plane): Plane ABH.
        BCH (Plane): Plane BCH.
        CDJ (Plane): Plane CDJ.
        DEK (Plane): Plane DEK.
        EFL (Plane): Plane EFL.
        FAG (Plane): Plane FAG.
        iterable_edges (list): List of all edges of the prism.
        iterable_planes (list): List of all planes of the prism.
        iterable_facets (list): List of all facets of the prism.

    Methods:
        is_point_in_truncated_triangle(point): Checks if a point is inside the truncated triangle.
        intersection_with_plane(plane): Finds the intersection points of the prism with a given plane.
        intersection_with_basal(plane): Finds the intersection points of the prism with the basal planes.

    """
    def __init__(self, r, k, m=1):
        """
        Initializes a WCCoPrism object with the given parameters.

        Args:
        r (float): Length of the first segment.
        k (float): Length of the second segment.
        m (float, optional): Length of the third segment. Default is 1. Default a+2b=1.
        """
        
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
                    for plane in self.iterable_planes:
                        if (p in plane) and (point in plane):
                            polygon.append(p)
                            del intersections[i]
                            point = p
                            flag = True
                            break
                    if flag:
                        break
                if uu_kostyl > 100:
                    #print('i had a problem with precision')
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
    '''
    Generates a path based on the parameters r, k, and n.

    Args:
        r (float): The parameter r.
        k (float): The parameter k.
        n (int, optional): The number of steps. Default is 100000.

    Returns:
        str: A string representing the path based on the input parameters.
    '''
    return "r=" + str(r) + "_k=" + str(k) + "_n=" + str(n)


def secant_polygon(polyhedron, way, params=None):
    '''
    Generates a plane inside the polyhedron using the specified method and returns a list of vertices of the secant polygon.

    Args:
        polyhedron: An object representing the polyhedron.
        way: A method for generating the plane.
        params (dict, optional): Parameters for the method, if required. Default is None.

    Returns:
        list: A list of vertices of the secant polygon.

    Notes:
        It is recommended to choose a way that fully tiles the polyhedron.
        The polyhedron object must have the method 'intersection_with_plane' implemented.
    '''
    if params:
        plane = way(**params)
    else:
        plane = way()
    # Select only planes that have a non-empty intersection with the polyhedron
    points_intersection = polyhedron.intersection_with_plane(plane)
    return points_intersection
        

def norm_dict(char_dict, step):
    '''
    Normalizes a characteristic provided as a dictionary to obtain the probability density from a histogram.

    Parameters:
        char_dict: dict
            Dictionary representing the characteristic histogram.
        step: float
            Step size used in the histogram.

    Returns:
        dict
            Normalized dictionary representing the probability density.
    '''
    norm_const = step*sum(char_dict.values())
    return {key: value/norm_const for key, value in char_dict.items()}


def getMinVolEllipse(P=None, tolerance=0.0001):
    """ Find the minimum volume ellipsoid which holds all the points
    code from https://github.com/minillinim/ellipsoid/blob/master/ellipsoid.py

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
    Computes the rotation matrix from one plane containing vectors Ax and Ay to another plane containing vectors Bx and By.

    Parameters:
        Ax: numpy.ndarray
            Vector representing the x-axis of the first plane.
        Ay: numpy.ndarray
            Vector representing the y-axis of the first plane.
        Bx: numpy.ndarray
            Vector representing the x-axis of the second plane.
        By: numpy.ndarray
            Vector representing the y-axis of the second plane.

    Returns:
        numpy.ndarray
            Rotation matrix from the first plane to the second plane.

    Notes:
        yes!! yes!! science!!
    '''
    #make Ax атв Ay identity and orthogonal
    Ax = Ax/np.linalg.norm(Ax)
    Az = (np.cross(Ax, Ay))/np.linalg.norm(np.cross(Ax, Ay))
    Ay = (np.cross(Ax, Az))/np.linalg.norm(np.cross(Ax, Az))
    
    Bx = Bx/np.linalg.norm(Bx)
    Bz = (np.cross(Bx, By))/np.linalg.norm(np.cross(Bx, By))
    By = (np.cross(Bx, Bz))/np.linalg.norm(np.cross(Bx, Bz))
    
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
    Rotates a set of 3D points to align with the Oxy plane. points should be in the same plane

    Parameters:
        points: numpy.ndarray
            Array of shape (N, 3) containing the 3D coordinates of points.

    Returns:
        numpy.ndarray
            Array of shape (N, 2) containing the rotated 2D coordinates of points aligned with the Oxy plane.
        numpy.ndarray
            Rotation matrix used for the transformation.
        float
            z-coordinate of the rotated points.
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
    Rotates a set of 2D points back to 3D coordinates.

    Parameters:
        rotated_polygon: numpy.ndarray
            Array of shape (N, 2) containing the rotated 2D coordinates of points.
        rotation_matrix: numpy.ndarray
            Rotation matrix used for the transformation.
        z: float
            z-coordinate of the points after rotation.

    Returns:
        numpy.ndarray
            Array of shape (N, 3) containing the rotated 3D coordinates of points.
    '''
    polygon = np.zeros((len(rotated_polygon), 3))
    polygon[:, :2] = rotated_polygon
    polygon[:, 2] = z
    polygon = np.array([np.dot(rotation_matrix.T, point) for point in polygon])
    return polygon


def is_convex_polygon_2D(points):
    """
    Checks whether a 2D polygon is convex.
    
    Arguments:
    points (list of tuples): List of points representing the vertices of the polygon. Each point is a tuple (x, y).
    
    Returns:
    bool: True if the polygon is convex, False otherwise.
    
    Note:
    A polygon must have at least 3 vertices to be considered convex.
    """

    n = len(points)
    if n < 3:
        return False  # A polygon must contain at least 3 vertices

    # Check the orientation for each triplet of consecutive vertices
    # If the sign of the angle between consecutive edges changes, the polygon is not convex
    def cross_product(p1, p2, p3):
        return (p2[0] - p1[0]) * (p3[1] - p2[1]) - (p2[1] - p1[1]) * (p3[0] - p2[0])

    direction = None
    for i in range(n):
        p1 = points[i]
        p2 = points[(i + 1) % n]
        p3 = points[(i + 2) % n]
        cp = cross_product(p1, p2, p3)
        if cp == 0:
            continue  # Points are collinear, continue checking
        if direction is None:
            direction = cp > 0
        elif direction != (cp > 0):
            return False  # Sign of the angle changed, polygon is not convex
    return True 


def semiaxes_from_polygon_3D(polygon):
    """
    Calculates the lengths of the semiaxes of the minimum-volume ellipsoid enclosing a 3D polygon.

    Arguments:
    polygon (list of tuples): List of points representing the vertices of the 3D polygon.

    Returns:
    tuple: A tuple containing the lengths of the semiaxes of the minimum-volume ellipsoid. 
    The first value represents the length of the minor semiaxis, and the second value represents the length of the major semiaxis.

    Note:
    The polygon should contain at least 3 non-collinear points.
    """
    polygonOxy, _, _ = rotate_points_to_Oxy(polygon)
    _, radii, _ = getMinVolEllipse(polygonOxy)
    return min(radii), max(radii)


def perimeter_3D(points):
    """
    Calculates the perimeter of a 3D polygon defined by its vertices.

    Arguments:
    points (list of tuples or numpy array): List of points representing the vertices of the 3D polygon.

    Returns:
    float: The perimeter of the 3D polygon.
    """
    polygon = np.array([np.array(list(i)) for i in points])
    per = 0.
    n = len(points)
    for i in range(n):
        per += np.linalg.norm(polygon[(i+1) % n]-polygon[i])
    return per


def minmax_feret(polygon):
    """
    Calculates the minimum and maximum Feret diameters of a 3D polygon.

    Arguments:
    polygon (numpy array): Array containing the vertices of the 3D polygon.

    Returns:
    tuple: A tuple containing the minimum and maximum Feret diameters of the 3D polygon.
    """
    _, ret_ix = antipodal_pairs(polygon)

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
    """
    Calculates various characteristics of a 3D section obtained by intersecting a polyhedron.

    Arguments:
    secant_polygon (list of tuples or numpy array): List of points representing the vertices of the 3D section.

    Returns:
    tuple: A tuple containing the following characteristics:
        - angles (list of float): List of angles in the section.
        - area (float): Area of the section.
        - perimeter (float): Perimeter of the section.
        - a_semi (float): Semi-major axis of the ellipse fitted to the section.
        - b_semi (float): Semi-minor axis of the ellipse fitted to the section.
        - min_feret (float): Minimum Feret diameter of the section.
        - max_feret (float): Maximum Feret diameter of the section.
        - aspect_ratio (float): Aspect ratio of the fitted ellipse (min_feret / max_feret).
        - sphericity (float): Sphericity of the section.
    """

    if len(secant_polygon) > 2:
        angles = get_angles(secant_polygon)
        area = area_of_polygon_3d(secant_polygon)
        perim = perimeter_3D(secant_polygon)
        a_semi, b_semi = semiaxes_from_polygon_3D(secant_polygon)
        min_feret, max_feret = minmax_feret(secant_polygon)
        aspect_ratio = min_feret/max_feret
        sphericity = 2*(np.pi*area)**0.5/perim
        
    else:
        angles = []
        area = 0.00
        a_semi, b_semi = 0.0, 0.0
        perim = 0.,
        min_feret, max_feret = 0., 0.
        aspect_ratio = 0.
        sphericity = 0.
  
    return angles, area, perim, a_semi, b_semi, min_feret, max_feret, aspect_ratio, sphericity


def calculate_all_distributions(polyhedron, way, PATH, type_figure='parallelepiped', n=10000, params=None):
        """
        Generate n planes using the method way and return the density distribution
        of geometric characteristics of the sections: angles, areas, perimeters, a_semi, b_semi, min_feret, max_feret, aspect ratio, and sphericity.
        
        If the directory 'densities' does not exist in the PATH, it is created, and the computed distributions are saved there.
        
        Parameters:
        polyhedron : Polyhedron
            An instance of a polyhedron class representing the shape in which to generate the sections.
        way : str
            The method to generate the planes. The specifics of this parameter depend on the implementation details.
        PATH : str
            File path to save the intercept lengths(computed distributions are saved in the folder 'densities' in PATH).
        type_figure : str, optional
            Type of figure to generate ('parallelepiped'/'triangular prism'/'hex prism').
        n : int, optional
            Number of sections to generate. Default is 10000.
        params : dict
            Parameters of the polyhedron.
            
        Returns:
            angles, areas, perimeters, a_semi, b_semi, min_feret, max_feret, aspect ratio, and sphericity.
        """

        PATH = os.path.join(PATH, f"densities")
        os.makedirs(PATH, exist_ok=True)


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

        elif type_figure == 'triangular prism':

            if params:
                side = params['side']
            else:
                side = 1.
                params = {'side': 1.}

            path_ = os.path.join(PATH, f"tri-prism_side={side}")

            angles = [0]*181
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
        
        
        elif type_figure == 'hex prism':

            r = params['r']
            k = params['k']

            path_ = os.path.join(PATH, f"hex-prism_r={r}k={k}")

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
        else:
            print('type figure should be "parallelepiped", "triangular prism" or "hex prism"')
            return None


        num_errors = 0

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
    Generates a random point uniformly within the triangle defined by the vertices A1, A2, and A3.
    
    Parameters:
    A1, A2, A3: array-like
        2D or 3D coordinates of the triangle vertices.

    Returns:
    Point
        A point generated inside the triangle with coordinates.
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
    Check if point P is inside the polygon.
    
    Parameters:
    P: array-like or np.ndarray
        Coordinates of the point to check.
    polygon: np.ndarray
        Coordinates of the polygon vertices.
    
    Returns:
    bool
        True if point is inside the polygon, False otherwise.
    '''
    S = area_of_polygon_3d(polygon) #area of polygon
    S_ = 0 #area of triangles PAB, where A, B - neighbour points of polygon
    
    N_edges = len(polygon)
    for i in range(N_edges-1):
        S_ += area_of_triangle_3d([P, polygon[i], polygon[i+1]])
        
    S_ += area_of_triangle_3d([P, polygon[N_edges-1], polygon[0]])
    
    return abs(S-S_) < 0.0000001


def proportion_S(triangles):
    '''
    Calculate the proportion of each triangle's area over the sum of all triangle areas.
    
    Parameters:
    triangles: np.ndarray
        An array of triangles where each triangle is represented by its vertices.
    
    Returns:
    np.ndarray
        An array of proportions representing the area of each triangle over the sum of all triangle areas.
    '''
    n_triangles = triangles.shape[0]
    proportions = np.zeros(n_triangles)
    
    for j in range(n_triangles):
        proportions[j] = area_of_triangle_3d(triangles[j])

    return proportions/proportions.sum()


def generations_point_in_polygon(vectors):
    '''
    Generate a random point inside a triangle defined by points A1, A2, and A3.
    
    Parameters:
    A1, A2, A3: np.ndarray or list
        2D or 3D coordinates of the triangle vertices.
    
    Returns:
    np.ndarray
        A random point inside the triangle.
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
    Generates a random plane and projects the vertices of the cube onto this plane.
    
    Parameters:
    cube : Cube (default is a unit cube)
        An instance of the Cube class representing the cube to be projected.
        
    Returns:
    tuple : (projected_points, hull)
        projected_points : np.ndarray
            Array of points projected onto the plane.
        hull : scipy.spatial.ConvexHull
            Convex hull of the projected points.
    '''
    # Generate a random normal vector
    fi = 2*np.pi*uniform(0, 1)
    teta = np.arccos(2*uniform(0, 1) - 1)
    n = Vector(np.sin(teta)*np.cos(fi), 
               np.sin(teta)*np.sin(fi),
               np.cos(teta))
    
    # Define a plane with the normal vector passing through the origin
    plane = Plane(Point(0., 0., 0.), n)
    
    # Project the cube vertices onto the plane
    polygon_points_not_ordered = []
    for M in cube.iterable_points:
        line_MP = Line(M, n)
        # P - projection of M on the plane
        P = intersection(line_MP, plane)
        polygon_points_not_ordered.append([P[0], P[1], P[2]])
    
    polygon_points_not_ordered = np.array(polygon_points_not_ordered)
    hull = ConvexHull(polygon_points_not_ordered[:, :2])
    
    return polygon_points_not_ordered, hull


def generate_linear_intersept(polyhedron):
    '''
    Generates the points and a length of a random linear intercept in a polyhedron.
    
    Parameters:
    polyhedron : Polyhedron
        An instance of a polyhedron class representing the shape in which to generate the intercept.
        
    Returns:
    tuple : (endpoints, length)
        endpoints : tuple of Points
            The two endpoints of the linear intercept.
        length : float
            The length of the linear intercept.
    '''
    # Generate a random normal vector
    fi = 2 * np.pi * uniform(0, 1)
    teta = np.arccos(2 * uniform(0, 1) - 1)
    n = Vector(np.sin(teta) * np.cos(fi), 
               np.sin(teta) * np.sin(fi),
               np.cos(teta))
    
    # Define a plane with the normal vector passing through the origin
    plane = Plane(Point(0., 0., 0.), n)
    
    # Project the polyhedron vertices onto the plane
    polygon_points_not_ordered = []
    for M in polyhedron.iterable_points:
        line_MP = Line(M, n)
        # P - projection of M on the plane
        P = intersection(line_MP, plane)
        polygon_points_not_ordered.append([P[0], P[1], P[2]])
    
    polygon_points_not_ordered = np.array(polygon_points_not_ordered) 
    hull = ConvexHull(polygon_points_not_ordered[:, :2]) # Compute convex hull of the projection
    polygon = polygon_points_not_ordered[hull.vertices]  # Extract vertices of the convex hull
    
    Pns = generations_point_in_polygon(polygon) # Generate a point in the projected polygon
    linear_intersept_line = Line(Pns, n) # Line through Pns in the direction of n
    
    two_points = [] # List to store intersection points with polyhedron facets
    
    for i, pl in enumerate(polyhedron.iterable_planes):
        # Check the intersection of generated line and the polyhedron's facets
        X = intersection(linear_intersept_line, pl)
        # Check if X is inside the facet
        if is_point_in_polygon(X, from_list_points_to_array(polyhedron.iterable_facets[i])):
            two_points.append(X)
    if len(two_points) != 2:
        #print('Something went wrong, not 2 points')
        return None, None
    P1 = two_points[0]
    P2 = two_points[1]
    lin = P1 - P2
    
    return (P1, P2), (lin[0]**2 + lin[1]**2 + lin[2]**2)**0.5


def calculate_linear_intercept(polyhedron, PATH, type_figure='parallelepiped', n=10000, params=None):
    """
    Calculate and save lengths of random linear intercepts in a polyhedron.
    If the directory 'densities' does not exist in the PATH, it is created, and the computed distributions are saved there.
    
    Parameters:
    polyhedron : Polyhedron
        An instance of a polyhedron class representing the shape in which to generate the intercepts.
    PATH : str
        File path to save the intercept lengths(computed distributions are saved in the folder 'densities' in PATH).
    type_figure : str
        Type of figure to generate ('parallelepiped'/'triangular prism'/'hex prism').
    n : int
        Number of intercepts to generate. Default is 10000.
    params : dict
        Parameters of the polyhedron.
        
    Returns:
        intercept
    """
    #PATH = os.path.join(os.getcwd(), 'densities')

    PATH = os.path.join(PATH, f"densities")
    os.makedirs(PATH, exist_ok=True)

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

        path_ = os.path.join(PATH, f"parallelepiped_a={a_side}b={b_side}c={c_side}")

        intercept_dict = {} #пустой словать с ключами от 0.00 до a^2+b^2+c^2 #можно меньше
        step = 0.01
        current_key = 0.00
        while current_key <= 10*(a_side**2 + b_side**2 + c_side**2):
            intercept_dict.setdefault(round(current_key, 2), 0)
            current_key += step

    elif type_figure == 'triangular prism':

        if params:
            side = params['side']
        else:
            side = 1.
            params = {'side': 1.}

        #path_ = PATH + r'\tri-prism_side=' + str(side)
        path_ = os.path.join(PATH, f"tri-prism_side={side}")

        intercept_dict = {} #пустой словать с ключами от 0.00 до a^2+b^2+c^2 #можно меньше
        step = 0.01
        current_key = 0.00
        while current_key <= 10*(side*side):
            intercept_dict.setdefault(round(current_key, 2), 0)
            current_key += step

    elif type_figure == 'hex prism':

        r = params['r']
        k = params['k']

        path_ = os.path.join(PATH, f"hex-prism_r={r}k={k}")

        intercept_dict = {} #пустой словать с ключами от 0.00 до a^2+b^2+c^2 #можно меньше
        step = 0.01
        current_key = 0.00
        while current_key <= 10.:
            intercept_dict.setdefault(round(current_key, 2), 0)
            current_key += step
    else:
        print('type figure should be "parallelepiped", "triangular prism" or "hex prism"')
        return None

    try:
        intercept = np.load(os.path.join(path_, f'intercept_N={n}.npy'), allow_pickle=True).item()
    except:
        num_errors = 0   
        for i in range(n):
            _, temp = generate_linear_intersept(polyhedron=polyhedron)
            if temp is not None:
                intercept_dict[round(temp, 2)] += 1
            else:
                num_errors += 1

        print('precision is ', 1-num_errors/n)

        intercept = norm_dict(intercept_dict, step=0.01)

        os.makedirs(path_, exist_ok=True)
        np.save(os.path.join(path_, f'intercept_N={n}.npy'), intercept)

    return intercept