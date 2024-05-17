from matplotlib import pyplot as plt
import numpy as np
from .funcs import secant_polygon
from .funcs import Parallelepiped, RegularTriangularPrism, WCCoPrism
from .funcs import make_random_plane_into_tetra, make_random_plane_into_prism, make_random_plane_into_wcco, generate_linear_intersept

from .funcs import area_of_polygon_3d, from_list_points_to_array, characteristic_of_section

def plot_secant_polygon(plane, type_figure='cube', params=None):
    '''
    Plots a 3D graph of a specified type_figure and its random section.

    Parameters:
        plane: Plane
            A plane that intersect polyhedron
        type_figure: str, optional
            Type of the figure. Default is 'cube'.
        params: dict, optional
            Parameters needed to define the polyhedron.

    Returns:
        None

    Notes:
       types figures and params should be:
       'cube': 'side',
       'parallelepiped': 'side_a', 'side_b', 'side_c',
       'triangular prism': 'side',
       'hex prism': 'r', 'k'
    '''

    if type_figure == 'cube':
        params = {'side_a': params['side'],
                  'side_b': params['side'],
                  'side_c': params['side']}
        
        type_figure = 'parallelepiped'

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if type_figure=='parallelepiped':
        polyhedron = Parallelepiped(**params)
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)
        ]
        secant = secant_polygon(polyhedron, make_random_plane_into_tetra, params)
        plt.title('Parallelepiped with a={}, b={}, c={}'.format(params['side_a'], params['side_b'], params['side_c']))


    if type_figure=='triangular prism':
        polyhedron = RegularTriangularPrism(**params)
        edges = [
            (0, 1), (1, 2), (2, 0),
            (3, 4), (4, 5), (5, 3),
            (0, 3), (1, 4), (2, 5)
        ]
        secant = secant_polygon(polyhedron, make_random_plane_into_prism, params)
        plt.title('Regular triangular prism with side={}'.format(params['side']))


    if type_figure=='hex prism':
        polyhedron = WCCoPrism(**params)
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
            (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6),
            (0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)
        ]
        secant = secant_polygon(polyhedron, make_random_plane_into_wcco, params)
        plt.title('Hex prism with r={}, k={}'.format(params['r'], params['k']))
        
    vertices = from_list_points_to_array(polyhedron.iterable_points)
    # edges of polyhedron
    for edge in edges:
        x = [vertices[edge[0]][0], vertices[edge[1]][0]]
        y = [vertices[edge[0]][1], vertices[edge[1]][1]]
        z = [vertices[edge[0]][2], vertices[edge[1]][2]]
        ax.plot(x, y, z, color='purple')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

    # Calculate characteristics of the section
    angles, area, perim, a_semi, b_semi, min_feret, max_feret, aspect_ratio, sphericity = characteristic_of_section(secant)

    # Print the results
    print('Characteristics of the section:')
    print("Angles:", angles)
    print("Area:", area)
    print("Perimeter:", perim)
    print("Semi-major axis:", a_semi)
    print("Semi-minor axis:", b_semi)
    print("Min Feret diameter:", min_feret)
    print("Max Feret diameter:", max_feret)
    print("Aspect ratio:", aspect_ratio)
    print("Sphericity:", sphericity)        
    

    # lines between vertices of secant polygon
    for i in range(len(secant)):
        x1, y1, z1 = secant[i]
        x2, y2, z2 = secant[(i + 1) % len(secant)]
        ax.plot([x1, x2], [y1, y2], [z1, z2], color='green')

    # add vertices of secant polygon
    for vertex in secant:
        x, y, z = vertex
        ax.scatter(x, y, z, color='blue')
    plt.show()


def plot_intercept(type_figure='cube', params=None):
    '''
    Plots a 3D graph of a specified type_figure and its random section.

    Parameters:
        plane: Plane
            A plane that intersect polyhedron
        type_figure: str, optional
            Type of the figure. Default is 'cube'.
        params: dict, optional
            Parameters needed to define the polyhedron.

    Returns:
        None

    Notes:
       types figures and params should be:
       'cube': 'side',
       'parallelepiped': 'side_a', 'side_b', 'side_c',
       'triangular prism': 'side',
       'hex prism': 'r', 'k'
    '''

    if type_figure == 'cube':
        params = {'side_a': params['side'],
                  'side_b': params['side'],
                  'side_c': params['side']}
        
        type_figure = 'parallelepiped'

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if type_figure=='parallelepiped':
        polyhedron = Parallelepiped(**params)
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7)
        ]
        
        plt.title('Parallelepiped with a={}, b={}, c={}'.format(params['side_a'], params['side_b'], params['side_c']))


    if type_figure=='triangular prism':
        polyhedron = RegularTriangularPrism(**params)
        edges = [
            (0, 1), (1, 2), (2, 0),
            (3, 4), (4, 5), (5, 3),
            (0, 3), (1, 4), (2, 5)
        ]
        plt.title('Regular triangular prism with side={}'.format(params['side']))


    if type_figure=='hex prism':
        polyhedron = WCCoPrism(**params)
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
            (6, 7), (7, 8), (8, 9), (9, 10), (10, 11), (11, 6),
            (0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)
        ]
        plt.title('Hex prism with r={}, k={}'.format(params['r'], params['k']))
        

    secant, length  = generate_linear_intersept(polyhedron)
    vertices = from_list_points_to_array(polyhedron.iterable_points)
    # edges of polyhedron
    for edge in edges:
        x = [vertices[edge[0]][0], vertices[edge[1]][0]]
        y = [vertices[edge[0]][1], vertices[edge[1]][1]]
        z = [vertices[edge[0]][2], vertices[edge[1]][2]]
        ax.plot(x, y, z, color='purple')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

    print("Length of linear intercept :", length)

    # intercept
    x1, y1, z1 = secant[0]
    x2, y2, z2 = secant[1]
    ax.plot([x1, x2], [y1, y2], [z1, z2], color='green')

    # add vertices of secant polygon
    for vertex in secant:
        x, y, z = vertex
        ax.scatter(x, y, z, color='blue')
    plt.show()


def ellipse_point(center, radii, rotation_matrix, num_points=100):
    '''
    Returns points of an ellipse with specified parameters.

    Parameters:
        center: tuple
            The coordinates of the center of the ellipse.
        radii: tuple
            The lengths of the semi-axes of the ellipse.
        rotation_matrix: numpy array
            The rotation matrix used to rotate the ellipse.
        num_points: int, optional
            The number of points to generate. Default is 100.

    Returns:
        numpy array
            An array of (x, y) coordinates representing points on the ellipse.
    '''
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
    
    return np.hstack([np.expand_dims(x_rotated, 0).T, np.expand_dims(y_rotated, 0).T])



def draw_ellipse(center, radii, rotation_matrix, num_points=100):
    '''
    Plots an ellipse with specified parameters.

    Parameters:
        center: tuple
            The coordinates of the center of the ellipse.
        radii: tuple
            The lengths of the semi-axes of the ellipse.
        rotation_matrix: numpy array
            The rotation matrix used to rotate the ellipse.
        num_points: int, optional
            The number of points to generate for the ellipse. Default is 100.
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

    