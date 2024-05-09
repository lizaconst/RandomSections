from matplotlib import pyplot as plt
import numpy as np


def ellipse_point(center, radii, rotation_matrix, num_points=100):
    '''
    возвращает точки эллипса с заданными параметрами
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

    