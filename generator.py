import random
import numpy as np
from math import log, sqrt

def normaly_dispersed_number_generator(ExpectedValue, stdDev, n):
    a = np.zeros(n)
    for i in range(n):
        a[i] = generate_normaly_dispersed_number(ExpectedValue, stdDev)
    return a

import random
from math import sqrt, log

def generate_normaly_dispersed_number(ExpectedValue, stdDev):
    """
    Generates a normally dispersed random number using the Box-Muller transform.

    Parameters:
    - ExpectedValue (float): The expected value (mean) of the normal distribution.
    - stdDev (float): The standard deviation of the normal distribution.

    Returns:
    - float: A random number generated from the normal distribution.

    """
    while True:
        u = 2 * random.random() - 1
        v = 2 * random.random() - 1
        s = u ** 2 + v ** 2
        if 0 < s < 1:
            break
    r = (sqrt(-2 * log(s) / s))
    return ExpectedValue + stdDev * random.choice([u, v]) * r

def generate_input_matrix(rows, cols, lower, upper):
    return np.random.uniform(lower, upper, (rows, cols))

def generate_input_vector(n, lower, upper, zeros):
    a = np.random.randint(lower, upper, n)
    zeros_count = np.count_nonzero(a == 0)
    while zeros_count < zeros:
        i = random.randint(0, n - 1)
        if a[i] != 0:
            a[i] = 0
            zeros_count += 1
    return a

def fetalon(indeces, row, operations):
    return sum([operations[i](row) * indeces[i] for i in range(len(operations))])

import numpy as np

def generate_input_set(input_size, repeats_per_exp, d, indeces, mat, exp_split, operations):
    """
    Generate an input set for experiments.

    Args:
        input_size (int): The size of the input set.
        repeats_per_exp (int): The number of repeats per experiment.
        d (float): The dispersion value.
        indeces (list): The list of indices.
        mat (numpy.ndarray): The matrix.
        exp_split (int): The split index for the output.
        operations (list): The list of operations.

    Returns:
        numpy.ndarray: The generated input set.
    """
    y0 = np.zeros((input_size, repeats_per_exp))
    for i in range(repeats_per_exp):
        for j in range(input_size):
            y0[j, i] = fetalon(indeces, mat[j], operations)

    Error = [normaly_dispersed_number_generator(0, d, repeats_per_exp) for i in range(input_size)]
    y = y0 + Error

    return y[:, :exp_split], y[:, exp_split:]
