import numpy as np
import os


def read_file_into_lines(input_file_path):
    with open(input_file_path) as fp_read:
        lines = []
        for line in fp_read:
            line = line.rstrip('\n')
            line = line.strip()
            lines.append(line)
    return lines


def sign(x):
    if x < 0:
        return -1
    elif x == 0:
        return 0
    else:
        return 1


def reverse_dict(input_dict):
    """
    requires mapping to be bijective
    :param input_dict:
    :return:
    """
    output_dict = {}
    for position, name in input_dict.items():
        output_dict[name] = position
    return output_dict


def generate_parabola(start_x, end_x, peak_x, peak_y, num_points=100):
    # Solve for the coefficients of the parabola
    A = np.array([[start_x ** 2, start_x, 1],
                  [peak_x ** 2, peak_x, 1],
                  [end_x ** 2, end_x, 1]])
    b = np.array([0.5, peak_y, 0.5])
    coef = np.linalg.solve(A, b)

    # Generate x values
    x = np.linspace(start_x, end_x, num_points)
    # Calculate y values for the parabola
    y = coef[0] * x ** 2 + coef[1] * x + coef[2]
    # Return x, y as tuples with precision to the sixth decimal point
    return np.array([(round(x_val, 6), round(y_val, 6)) for x_val, y_val in zip(x, y)])


def generate_circle(start_x, peak_y, circle_size_multiplier, num_points=100):
    diameter = abs(peak_y - 0.5) * circle_size_multiplier
    radius = diameter / 2
    is_upward_direction = (peak_y > 0.5)

    # Calculate the center of the circle
    center_x = start_x
    if is_upward_direction:
        center_y = 0.5 + radius  # Fixed y-coordinate
        # Generate angles for the circle (starting from the bottom and going clockwise)
        angles = np.linspace(3 * np.pi / 2, -np.pi / 2, num_points)
    else:
        center_y = 0.5 - radius
        # Generate angles for the circle (starting from the top and going counterclockwise)
        angles = np.linspace(np.pi / 2, -3 * np.pi / 2, num_points)

    # Generate x and y coordinates for the circle
    x = center_x + radius * np.cos(angles)
    y = center_y + radius * np.sin(angles)

    # Combine x and y coordinates into points
    points = np.array([[round(x_val, 6), round(y_val, 6)] for x_val, y_val in zip(x, y)])

    # Make sure the first and last points are equal to the starting point
    if points[0][0] != start_x or points[-1][0] != start_x:
        raise RuntimeError('circle coordinate formation error')

    return points


def get_metadata_file_path(filename):
    karutils_dir = os.path.dirname(__file__)
    metadata_dir = os.path.join(karutils_dir, 'Metadata')
    return os.path.join(metadata_dir, filename)


def round_half_up(n):
    import math
    return int(math.floor(n + 0.5))


def convert_chrom(i_chrom):
    if i_chrom.upper() == 'X':
        return 23
    elif i_chrom.upper() == 'Y':
        return 24
    else:
        return int(i_chrom)
