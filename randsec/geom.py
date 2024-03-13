from decimal import Decimal
from fractions import Fraction

import math
import numpy as np

def unify_types(items):
    """Promote all items to the same type. The resulting type is the
    "most valueable" that an item already has as defined by the list
    (top = least valueable):
    - int
    - float
    - decimal.Decimal
    - fractions.Fraction
    - user defined
    """
    type_values = {
        Fraction: 1,
        Decimal: 2,
        float: 3,
        int: 4,
    }
    types = []
    for item in items:
        for type_, value in type_values.items():
            if isinstance(item, type_):
                types.append((value, type_))
                break
        else:
            types.append((0, type(item)))
    result_type = min(types)[1]
    return [result_type(i) for i in items]



#solver.py

def shape(m):
    if not m:
        return (0, 0)
    return (len(m), len(m[0]))

def null(f):
    return abs(f) < 1e-10

def nullrow(r):
    return all(map(null, r))

def find_pivot_row(m):
    candidates = []
    for i, row in enumerate(m):
        # Only rows where the pivot element is not zero can be used
        if row[0] != 0:
            candidates.append((abs(row[0]), i))
    if not candidates:
        return None
    # We use the one with the biggest absolute value
    return max(candidates)[1]

def gaussian_elimination(m):
    """Return the row echelon form of m by applying the gaussian
    elimination"""
    # Shape of the matrix
    M, N = shape(m)
    for j in range(N-1):
        # We ignore everything above the jth row and everything left of
        # the jth column (we assume they are 0 already)
        pivot = find_pivot_row([row[j:] for row in m[j:]])
        if pivot is None:
            continue
        # find_pivot_row returns the index relative to j, so we need to
        # calculate the absolute index
        pivot += j
        # Swap the rows
        m[j], m[pivot] = m[pivot], m[j]
        # Note that the pivot row is now m[j]!
        # Eliminate everything else
        for i in range(j + 1, M):
            factor = m[i][j] / m[j][j] * -1
            # Multiply the pivot row before adding them
            multiplied_row = [factor * x for x in m[j]]
            # Looks ugly, but we don't need numpy for it
            # Replace the ith row with the sum of the ith row and the
            # pivot row
            m[i] = [x + y for x, y in zip(m[i], multiplied_row)]
    # m shold now be in row echelon form
    return m

def solve(matrix):
    ref = gaussian_elimination(matrix)
    return Solution(ref)

def count(f, l):
    c = 0
    for i in l:
        if f(i):
            c += 1
    return c

def index(f, l):
    for i, v in enumerate(l):
        if f(v):
            return i
    raise ValueError("No item satisfies {}".format(f))

def first_nonzero(r):
    for i, v in enumerate(r):
        if not null(v):
            return i
    return len(r)

class Solution(object):
    """Holds a solution to a system of equations."""
    def __init__(self, s):
        self._s = s
        self.varcount = shape(s)[1] - 1
        # No solution, 0a + 0b + 0c + ... = 1 which can never be true
        self._solvable = not any(
            all(null(coeff) for coeff in row[:-1]) and not null(row[-1])
            for row in s
        )
        unique_equations = sum(1 for row in s if not nullrow(row))
        self.varargs = self.varcount - unique_equations
        self.exact =  self.varargs == 0

    def __bool__(self):
        return self._solvable
    __nonzero__ = __bool__

    def __call__(self, *v):
        if not self._solvable:
            raise ValueError("Has no solution")
        if len(v) != self.varargs:
            raise ValueError("Expected {} values, got {}".format(
                self.varargs, len(v)))
        v = list(v)
        vals = [None] * self.varcount
        # Scan for real solutions
        for i, row in enumerate(self._s):
            # Can't use .count here because we need null()
            # I miss Haskell lambdas :(
            if count(lambda i: not null(i), row[:-1]) == 1:
                # We can find a variable here
                var = index(lambda i: not null(i), row[:-1])
                vals[var] = row[-1] / row[var]
        # Fill in the rest with given values
        for i in reversed(range(len(vals))):
            if not v:
                break
            if vals[i] is None:
                vals[i] = v.pop()

        for i in reversed(range(len(self._s))):
            row = self._s[i]
            if nullrow(row):
                continue
            tbd = first_nonzero(row)
            s = sum(-1 * row[j] * vals[j] for j in range(tbd + 1, len(row) - 1))
            s += row[-1]
            vals[tbd] = s / row[tbd]
        return tuple(vals)
    


#body.py
    
class GeoBody(object):
    """A base class for geometric objects that provides some common
    methods to work with. In the end, everything is dispatched to
    sgl.calc.* anyway, but it sometimes feels nicer to write it like
    L1.intersection(L2)
    instead of
    intersection(L1, L2)
    """
    def intersection(self, other):
        #from calc import intersection
        return intersection(self, other)

    def distance(self, other):
        #from calc import distance
        return distance(self, other)

    def parallel(self, other):
        #from calc import parallel
        return parallel(self, other)

    def angle(self, other):
        #from calc import angle
        return angle(self, other)

    def orthogonal(self, other):
        #from calc import orthogonal
        return orthogonal(self, other)


#point.py
    
class Point(object):
    """Provides a basic Point in the 3D space"""
    @classmethod
    def origin(cls):
        """Returns the origin (0 | 0 | 0)"""
        return cls(0, 0, 0)
    o = origin

    def __init__(self, *args):
        """Point(a, b, c)
        Point([a, b, c]):
        The point with coordinates (a | b | c)

        Point(Vector):
        The point that you get when you move the origin by the given
        vector. If the vector has coordinates (a | b | c), the point
        will have the coordinates (a | b | c) (as easy as π).
        """
        if len(args) == 1:
            # Initialisation by Vector is also handled by this
            coords = args[0]
        elif len(args) == 3:
            coords = args
        else:
            raise TypeError("Point() takes one or three arguments, not {}"
                    .format(len(args)))
        self.x, self.y, self.z = unify_types(coords)

    def __repr__(self):
        return "Point({}, {}, {})".format(
                self.x,
                self.y,
                self.z,
                )

    def __hash__(self):
        return hash(("Point", self.x, self.y, self.z))

    def __add__(self, other):
        return Point(self.x + other.x, self.y + other.y, self.z + other.z)

    
    def __sub__(self, other):
        return Point(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, other):
        return Point(other*self.x, other*self.y, other*self.y)

    def __rmul__(self, other):
        return self * other

    def __eq__(self, other):
        """Checks if two Points are equal. Always use == and not 'is'!"""
        return (self.x == other.x and
                self.y == other.y and
                self.z == other.z)

    def __getitem__(self, item):
        return (self.x, self.y, self.z)[item]

    def __setitem__(self, item, value):
        setattr(self, "xyz"[item], value)

    def pv(self):
        """Return the position vector of the point."""
        #from vector import Vector
        return Vector(self.x, self.y, self.z)

    def moved(self, v):
        """Return the point that you get when you move self by vector v."""
        return Point(self.pv() + v)
    

#vector.py
class Vector(object):
    """Provides a basic vector"""
    @classmethod
    def zero(cls):
        """Returns the zero vector (0 | 0 | 0)"""
        return cls(0, 0, 0)

    def __init__(self, *args):
        """Vector(x, y, z)
        Vector([x, y, z]):
        A vector with coordinates (x | y | z)

        Vector(P1, P2):
        A vector going from point P1 to P2.
        """
        if len(args) == 3:
            # Initialising with 3 coordinates
            self._v = list(args)
        elif len(args) == 2:
            # Initialising from point A to point B
            A, B = args
            self._v = [
                B.x - A.x,
                B.y - A.y,
                B.z - A.z,
            ]
        elif len(args) == 1:
            # Initialising with an array of coordinates
            self._v = list(args[0])
        else:
            raise TypeError("Vector() takes one, two or three parameters, "
                            "not {}".format(len(args)))
        self._v = unify_types(self._v)

    def __hash__(self):
        return hash(("Vector",) + tuple(self))

    def __repr__(self):
        return "Vector({}, {}, {})".format(*self._v)
    
    def __eq__(self, other):
        return (self._v == other._v)

    def __add__(self, other):
        return Vector(x+y for x, y in zip(self, other))
    
    def __sub__(self, other):
        return Vector([x-y for x, y in zip(self, other)])

    def __mul__(self, other):
        if isinstance(other, Vector):
            return sum(x*y for x, y in zip(self, other))
        return Vector([x*other for x in self._v])

    def __rmul__(self, other):
        return self * other
    
    def __neg__(self):
        return self * -1

    def __getitem__(self, item):
        return self._v[item]

    def __setitem__(self, item, value):
        self._v[item] = value

    def cross(self, other):
        r"""Calculates the cross product of two vectors, defined as
        _   _   / x2y3 - x3y2 \
        x × y = | x3y1 - x1y3 |
                \ x1y2 - x2y1 /

        The cross product is orthogonal to both vectors and its length
        is the area of the parallelogram given by x and y.
        """
        a, b = self._v, other._v
        return Vector(
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0]
                )

    def length(self):
        """Returns |v|, the length of the vector."""
        return (self * self) ** 0.5
    __abs__ = length

    def parallel(self, other):
        """Returns true if both vectors are parallel."""
        #from solver import solve
        if self == Vector.zero() or other == Vector.zero():
            return False
        if self == other:
            return True
        # linear combination:
        # a * self + b * other = 0
        solution = solve([
            [self[0], other[0], 0],
            [self[1], other[1], 0],
            [self[2], other[2], 0],
        ])
        # Trivial solution is a = b = 0
        # if there are no other solutions, the vectors are not parallel!
        # otherwise there are infinitely many solutions and the vectors
        # are parallel.
        if solution.exact:
            return False
        return True

    def orthogonal(self, other):
        """Returns true if the two vectors are orthogonal"""
        return self * other == 0

    def angle(self, other):
        """Returns the angle (in radians) enclosed by both vectors."""
        try:
            return math.acos((self * other) / (self.length() * other.length()))
        except:
            print("i got a problem with arccos", self, other)
            if self.length() * other.length() == 0:
                return 0
            if (self * other) / (self.length() * other.length()) > 0:
                return 0
            else:
                return np.pi
       

    def normalized(self):
        """Return the normalized version of the vector, that is a vector
        pointing in the same direction but with length 1.
        """
        # Division is not defined, so we have to multiply by 1/|v|
        return float(1 / self.length()) * self
    unit = normalized


#line.py
    
class Line(GeoBody):
    """Provides a line in 3d space"""
    def __init__(self, a, b):
        """Line(Point, Point):
        A Line going through both given points.

        Line(Point, Vector):
        A Line going through the given point, in the direction pointed
        by the given Vector.

        Line(Vector, Vector):
        The same as Line(Point, Vector), but with instead of the point
        only the position vector of the point is given.
        """
        # We're storing the position vector, so if a point is given we
        # need to convert it first
        if isinstance(a, Point):
            a = a.pv()
        # Support vector
        self.sv = a
        if isinstance(b, Vector):
            self.dv = b
        elif isinstance(b, Point):
            # We just take the vector AB as the direction vector
            self.dv = b.pv() - self.sv

        if self.dv == Vector.zero():
            raise ValueError("Invalid Line, Vector(0 | 0 | 0)")

    def __repr__(self):
        return "Line({}, {})".format(self.sv, self.dv)

    def __contains__(self, point):
        """Checks if a point lies on a line"""
        v = point.pv() - self.sv
        return v.parallel(self.dv)

    def __eq__(self, other):
        """Checks if two lines are equal"""
        return Point(other.sv) in self and other.dv.parallel(self.dv)

    def parametric(self):
        """Returns (s, u) so that you can build the equation for the line
           _   _    _
        g: x = s + ru ; r e R
        """
        return (self.sv, self.dv)
    

#plane.py
    
class Plane(GeoBody):
    """A Plane (not the flying one)"""
    def __init__(self, *args):
        """Plane(Point, Point, Point):
        Initialise a plane going through the three given points.

        Plane(Point, Vector, Vector):
        Initialise a plane given by a point and two vectors lying on
        the plane.

        Plane(Point, Vector):
        Initialise a plane given by a point and a normal vector (point
        normal form)

        Plane(a, b, c, d):
        Initialise a plane given by the equation
        ax1 + bx2 + cx3 = d (general form).
        """
        if len(args) == 3:
            a, b, c = args
            if (isinstance(a, Point) and
                isinstance(b, Point) and
                isinstance(b, Point)):
                # for three points we just calculate the vectors AB
                # and AC and continue like we were given two vectors
                # instead
                vab = b.pv() - a.pv()
                vac = c.pv() - a.pv()
            elif (isinstance(a, Point) and
                  isinstance(b, Vector) and
                  isinstance(c, Vector)):
                vab, vac = b, c
            # We need a vector orthogonal to the two given ones so we
            # (the length doesn't matter) so we just use the cross
            # product
            vec = vab.cross(vac)
            self._init_pn(a, vec)
        elif len(args) == 2:
            self._init_pn(*args)
        elif len(args) == 4:
            self._init_gf(*args)
    
    def _init_pn(self, p, normale):
        """Initialise a plane given in the point normal form."""
        self.p = p
        self.n = normale

    def _init_gf(self, a, b, c, d):
        """Initialise a plane given in the general form."""
        # We need
        # 1) a normal vector -> given by (a, b, c)
        # 2) a point on the plane -> solve the equation and chose a
        #    "random" point
        solution = solve([[a, b, c, d]])
        self.n = Vector(a, b, c)
        self.p = Point(*solution(1, 1))

    def __eq__(self, other):
        """Checks if two planes are equal. Two planes can be equal even
        if the representation is different!
        """
        return self.p in other and self.parallel(other)

    def __contains__(self, other):
        """Checks if a Point lies on the Plane or a Line is a subset of
        the plane.
        """
        #from line import Line
        if isinstance(other, Point):
            #return other.pv() * self.n == self.p.pv() * self.n
            return abs(other.pv() * self.n - self.p.pv() * self.n) < 0.0001
        elif isinstance(other, Line):
            #return Point(other.sv) in self and self.parallel(other)
            return Point(other.sv) in self and other.dv*self.n < 0.0001

    def __repr__(self):
        return "Plane({}, {})".format(self.p, self.n)

    def point_normal(self):
        """Returns (p, n) so that you can build the equation
            _   _   
        E: (x - p) n = 0

        to describe the plane.
        """
        # That's the form we use to store the plane internally,
        # we don't have to calculate anything
        return (self.p.pv(), self.n)

    def general_form(self):
        """Returns (a, b, c, d) so that you can build the equation

        E: ax1 + bx2 + cx3 = d

        to describe the plane.
        """
        # Since this form is just the point-normal-form when you do the
        # multiplication, we don't have to calulate much here
        return (
            self.n[0],
            self.n[1],
            self.n[2],
            self.n * self.p.pv(),
        )

    def parametric(self):
        """Returns (u, v, w) so that you can build the equation
           _   _    _    _ 
        E: x = u + rv + sw ; (r, s) e R

        to describe the plane (a point and two vectors).
        """
        s = solve([list(self.n) + [0]])
        # Pick a first vector orthogonal to the normal vector
        # there are infinitely many solutions, varying in direction
        # and length, so just choose some values
        v = Vector(*s(1, 1))
        assert v.orthogonal(self.n)
        # Pick a second vector orthogonal to the normal vector and
        # orthogonal to the first vector (v)
        # again, there are infinitely many solutions, varying in length
        s = solve([
            list(self.n) + [0],
            list(v) + [0],
        ])
        w = Vector(*s(1))
        return (self.p.pv(), v, w)
    



#calc.py
    
def acute(rad):
    """If the given angle is >90° (pi/2), return the opposite angle"""
    if rad > 0.5 * math.pi:
        rad = math.pi - rad
    return rad

      
def intersection(a, b):
    """Return the intersection between two objects. This can either be
    - None (no intersection)
    - a Point (Line/Line or Plane/Line intersection)
    - a Line (Plane/Plane intersection)
    """
    if isinstance(a, Line) and isinstance(b, Line):
        # For the line-line intersection, we have to solve
        # s1 + λ u1 = t1 + μ v1
        # s2 + λ u2 = t2 + μ v2
        # s3 + λ u3 = t3 + μ v3
        # rearrange a bit, and you get
        solution = solve([
            [a.dv[0], -b.dv[0], b.sv[0] - a.sv[0]],
            [a.dv[1], -b.dv[1], b.sv[1] - a.sv[1]],
            [a.dv[2], -b.dv[2], b.sv[2] - a.sv[2]],
        ])
        # No intersection
        if not solution:
            return None
        # We get λ and μ, we need to pick one and plug it into the
        # right equation
        lmb, mu = solution()
        lmb = float(lmb)
        # could've chosen b.sv + mu * b.dv instead, it doesn't matter
        # as they will point (pun intended) to the same point.
        return Point(a.sv + lmb * a.dv)

    elif isinstance(a, Line) and isinstance(b, Plane):
        # the line can be contained in the plane, in this case the whole
        # line is the intersection
        if a in b:
            return a
        # if they are parallel, there is no intersection
        elif parallel(a, b):
            return None
        # Given the plane in general form, if we insert the line
        # coordinate by coordinate we get
        # a (s1 + μ u1) + b (s2 + μ u2) + c (s3 + μ u3) = d
        # where s is the support vector of the line
        #       u is the direction vector of the line
        #       μ is the parameter
        # rearrange and solve for the parameter:
        mu = (b.n * b.p.pv() - b.n * a.sv) / (b.n * a.dv)
        mu = float(mu)
        return Point(a.sv + mu * a.dv)
    elif isinstance(a, Plane) and isinstance(b, Line):
        return intersection(b, a)

    elif isinstance(a, Plane) and isinstance(b, Plane):
        # if you solve
        # a x1 + b x2 + c x3 = d
        # e x1 + f x2 + g x3 = h
        # you will get infinitely many solutions (if the planes are
        # intersecting). All those solutions are points on the 
        # intersection line. So we just chose two solutions, i.e.
        # two points, and lay a line through both of these.
        solution = solve([
            list(a.n) + [a.n * a.p.pv()],
            list(b.n) + [b.n * b.p.pv()],
        ])
        if not solution:
            return None
        # Choose two arbitrary points/solutions
        p1, p2 = Point(solution(1)), Point(solution(2))
        return Line(p1.pv(), p2.pv() - p1.pv())

    return NotImplemented


def parallel(a, b):
    """Checks if two objects are parallel. This can check
    - Line/Line
    - Plane/Line
    - Plane/Plane
    """
    if isinstance(a, Line) and isinstance(b, Line):
        return a.dv.parallel(b.dv)

    elif isinstance(a, Line) and isinstance(b, Plane):
        return a.dv.orthogonal(b.n)
    elif isinstance(a, Plane) and isinstance(b, Line):
        return parallel(b, a)
    
    elif isinstance(a, Plane) and isinstance(b, Plane):
        return a.n.parallel(b.n)

    return NotImplemented


def angle(a, b):
    """Returns the angle (in radians) between
    - Line/Line
    - Plane/Line
    - Plane/Plane
    """
    if isinstance(a, Line) and isinstance(b, Line):
        return acute(a.dv.angle(b.dv))

    elif isinstance(a, Line) and isinstance(b, Plane):
        rad = acute(a.dv.angle(b.n))
        # What we are actually calculating is the angle between
        # the normal of the plane and the line, but the normal
        # is 90° from the plane. So the actual angle between a plane
        # a line is 90° - that angle
        return 0.5 * math.pi - rad
    elif isinstance(a, Plane) and isinstance(b, Line):
        return angle(b, a)

    elif isinstance(a, Plane) and isinstance(b, Plane):
        return acute(a.n.angle(b.n))

    return NotImplemented


def orthogonal(a, b):
    """Checks if two objects are orthogonal. This can check
    - Line/Line
    - Plane/Line
    - Plane/Plane
    """
    if isinstance(a, Line) and isinstance(b, Line):
        return null(a.dv * b.dv)

    elif isinstance(a, Line) and isinstance(b, Plane):
        return a.dv.parallel(b.n)
    elif isinstance(a, Plane) and isinstance(b, Line):
        return orthogonal(b, a)

    elif isinstance(a, Plane) and isinstance(b, Plane):
        return a.n.orthogonal(b.n)
    
    return NotImplemented


def distance(a, b):
    """Returns the distance between two objects. This includes
    - Point/Point
    - Line/Point
    - Line/Line
    - Plane/Point
    - Plane/Line
    """
    if isinstance(a, Point) and isinstance(b, Point):
        # The distance between two Points A and B is just the length of
        # the vector AB
        return Vector(a, b).length()

    elif isinstance(a, Point) and isinstance(b, Line):
        # To get the distance between a point and a line, we place an
        # auxiliary plane P. P is orthogonal to the line and contains
        # the point. To achieve this, we just use the direction vector
        # of the line as the normal vector of the plane.
        aux_plane = Plane(a, b.dv)
        # We then calculate the intersection of the auxiliary plane and
        # the line
        foot = intersection(aux_plane, b)
        # And finally the distance between the point and the
        # intersection point, which can be reduced to a Point-Point
        # distance
        return distance(a, foot)
    elif isinstance(a, Line) and isinstance(b, Point):
        return distance(b, a)

    elif isinstance(a, Line) and isinstance(b, Line):
        # To get the distance between two lines, we just use the formula
        #        _   _    _
        # d = | (q - p) * n |
        # where n is a vector orthogonal to both lines and with length 1!
        # We can achieve this by using the normalized cross product
        normale = a.dv.cross(b.dv).normalized()
        return abs((b.sv - a.sv) * normale)

    elif isinstance(a, Point) and isinstance(b, Plane):
        # To get the distance between a point and a plane, we just take
        # a line that's orthogonal to the plane and goes through the
        # point
        aux_line = Line(a, b.n)
        # We then get the intersection point...
        foot = intersection(aux_line, b)
        # ...and finally the distance
        return distance(a, foot)
    elif isinstance(a, Plane) and isinstance(b, Point):
        return distance(b, a)

    elif isinstance(a, Line) and isinstance(b, Plane):
        if parallel(a, b):
            # If the line is parallel, every point has the same distance
            # to the plane, so we just pick one point and calculate its
            # distance
            return distance(Point(a.sv), b)
        # If they are not parallel, they will eventually intersect, so
        # the distance is 0
        return 0.0
    elif isinstance(a, Plane) and isinstance(b, Line):
        return distance(b, a)

    return NotImplemented