import numpy as np

def parse(n, filename):
    file = open(filename, "r")
    line = file.readline()
    row = [int(x) for x in line.split()]
    
    mat = np.zeros((n, len(row)-1))
    
    for i in range(n):
        line = file.readline()
        mat[i] = line.split()   
    return mat, row
