import numpy as np
import generator
import math 
from itertools import permutations

def least_sqares(mat, y):
    if mat[0].size == mat[:,0].size:
        if np.linalg.det(mat) == 0:
            raise Exception("Matrix is not invertible")
    ans = np.linalg.inv(np.transpose(mat) @ mat) @ np.transpose(mat) @ y
    return ans
    
def zsk(mat, y, teta):
    for i in range(y.size):
        t = 0
        for j in range(mat[0].size):
            sum = 0
            for k in range(mat[:,0].size):
                if mat[k][j] == 0:
                    return False


def sqrtOfSum(b): 
    n = len(b)
    sqrtOfSum = 0
    for i in range(n):
        sqrtOfSum += b[i]**2
    return math.sqrt(sqrtOfSum)

def comMeasure(b_og , b):
    n = len(b)
    sum = 0
    sumofb_og = sqrtOfSum(b_og)
    sumofb = sqrtOfSum(b)
    for i in range(n):
        sum += (b_og[i]/sumofb_og - b[i]/sumofb)**2
    return math.sqrt(sum)

def deploy_to_file(filename, data, d):
    f = open(filename, "a")
    f.write("Statistics:\n")
    f.write("Dispersion: " + str(d) + "\n")
    f.write("Average error: " + str(data[0]) + "\n")
    f.write("Max error: " + str(data[1]) + "\n")
    f.write("Min error: " + str(data[2]) + "\n")
    f.write("\n")
    f.close()

    
#input parameters

#polynom stucture
operations = [
    lambda  row: 1,
    lambda  row: row[0],
    lambda  row: row[1],
    lambda  row: row[2],
    lambda  row: row[3],
    lambda  row: row[0]*row[1],
    lambda  row: row[0]*row[3],
    lambda  row: row[1]**2,
    lambda  row: row[3]**3
    
]
    
    
#for input matrix

m_param = 4
lowerbound_m = 1
upperbound_m = 40
input_size = 100

#for input vector
inputvector_size = len(operations)
inputvector_lowerbound = -10
inputvector_upperbound = 10
inutvector_zeros = 3

#repeats
repeats_per_exp = 10
exp_count = 1000

#dispersion
dist_top = 10
d = [x for x in range(1, dist_top + 1)]

#exp split
exp_split = (input_size // 10) * 6

#end of input parameters
err = np.zeros((dist_top, exp_count))

for f in range(0, dist_top):

    for k in range(exp_count):
        mat = generator.generate_input_matrix(input_size, m_param, lowerbound_m, upperbound_m)
        indeces = generator.generate_input_vector(inputvector_size, inputvector_lowerbound, inputvector_upperbound, inutvector_zeros)
        y1, y2 = generator.generate_input_set(input_size, repeats_per_exp, d[f], indeces, mat, exp_split, operations)
        
        yA1 = np.array([sum(y1[i])/exp_split for i in range(input_size)])
        
        ml1 = np.array([[operations[j](mat[i]) for j in range(inputvector_size)] for i in range(input_size)])
        
        ans = least_sqares(ml1, yA1)    
        
        err[f,k] = comMeasure(indeces, ans)
        
             
        
    data = np.zeros(3)
    data[0] = np.mean(err[f])
    data[1] = np.max(err[f])
    data[2] = np.min(err[f])
    deploy_to_file("out.txt", data, d[f])

