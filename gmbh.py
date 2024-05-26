import numpy as np
import generator
import math 
from itertools import combinations

def least_sqares(mat, y):
    """
    Perform least squares regression to find the optimal solution.

    Args:
        mat (numpy.ndarray): The matrix of values.
        y (numpy.ndarray): The target values.

    Returns:
        numpy.ndarray: The optimal solution.
    """
    if mat[0].size == mat[:,0].size:
        if np.linalg.det(mat) == 0:
            raise Exception("Matrix is not invertible")
    ans = np.linalg.inv(np.transpose(mat) @ mat) @ np.transpose(mat) @ y
    return ans

def calcBestInd(mat, certanlyNonzeroIndex, comb, y, operations, yA):
    """
    Calculate the best combination of indices based on the given inputs.

    Args:
        mat (numpy.ndarray): The matrix of values.
        certanlyNonzeroIndex (list): The list of indices that are certainly non-zero.
        comb (list): The list of combinations of indices.
        y (numpy.ndarray): The target values.
        operations (list): The list of operations to be applied on the matrix.
        yA (numpy.ndarray): The matrix of transformed target values.

    Returns:
        list: The best combination of indices along with the certainly non-zero indices.
    """
    zskArr = [0]*len(comb)
    for cmbI in range(len(comb)):
        combinedIndex = certanlyNonzeroIndex + list(comb[cmbI])
        combinedIndex.sort()
        ml = np.array([[operations[j](mat[i]) for j in range(len(operations)) if j in combinedIndex] for i in range(len(mat))]) 
        tetaml = least_sqares(ml, yA)
        for k in range (len(y[0])):
            zskArr[cmbI] += zsk(ml, y[:, k], tetaml)
    min_zsk = min(zskArr)
    czsk = 0
    best_comb = None
    for i in range(len(zskArr)):
        if zskArr[i] <= 1.03 * min_zsk and zskArr[i] >= 0.97 * min_zsk:
            if best_comb is None or len(comb[i]) < len(best_comb):
                best_comb = comb[i]
                czsk = zskArr[i]
            elif len(comb[i]) == len(best_comb):
                if sum([abs(tetaml[j]) for j in comb[i]]) < sum([abs(tetaml[j]) for j in best_comb]):
                    best_comb = comb[i]
                    czsk = zskArr[i]
    return list(best_comb) + certanlyNonzeroIndex
    

def zsk(mat,  y, teta):
    """
    Calculate the ZSK (Squared Sum of Errors) value.

    Args:
        mat (numpy.ndarray): The matrix of values.
        y (numpy.ndarray): The target values.
        teta (numpy.ndarray): The solution vector.

    Returns:
        float: The ZSK value.
    """
    t = 0
    for i in range(len(y)):
        sumA = 0
        for j in range(len(mat[0])):
            sumA += mat[i][j] * teta[j]
        t += (y[i] - sumA) ** 2
    return t


def sqrtOfSum(b): 
    """
    Calculate the square root of the sum of squares of elements in a vector.

    Args:
        b (list): The input vector.

    Returns:
        float: The square root of the sum of squares.
    """
    n = len(b)
    sqrtOfSum = 0
    for i in range(n):
        sqrtOfSum += b[i]**2
    return math.sqrt(sqrtOfSum)

def comMeasure(b_og , b):
    """
    Calculate the coefficient of measurement.

    Args:
        b_og (list): The original vector.
        b (list): The calculated vector.

    Returns:
        tuple: The coefficient of measurement and the number of incorrect zeros.
    """
    zeroes_wrong = 0
    n = len(b)
    for i in range(n):
        if b[i] != 0 and b_og[i] == 0:
            zeroes_wrong += 1
    sum = 0
    sumofb_og = sqrtOfSum(b_og)
    sumofb = sqrtOfSum(b)
    for i in range(n):
        sum += (b_og[i]/sumofb_og - b[i]/sumofb)**2
    return math.sqrt(sum), zeroes_wrong

def deploy_to_file(filename, data, d):
    """
    Deploy the statistics to a file.

    Args:
        filename (str): The name of the file.
        data (list): The statistics data.
        d (int): The dispersion value.
    """
    f = open(filename, "a")
    f.write("Statistics:\n")
    f.write("Dispersion: " + str(d) + "\n")
    f.write("Average error: " + str(data[0]) + "\n")
    f.write("Max error: " + str(data[1]) + "\n")
    f.write("Min error: " + str(data[2]) + "\n")
    f.write("Average zeroes wrong: " + str(data[3]) + "\n")
    f.write("Max zeroes wrong: " + str(data[4]) + "\n")
    f.write("\n")
    f.close()

def find_val(val, arr):
    """
    Find the index of a value in a list.

    Args:
        val: The value to find.
        arr (list): The input list.

    Returns:
        int: The index of the value in the list, or -1 if not found.
    """
    for i in range(len(arr)):
        if arr[i] == val:
            return i 
    return -1

def find_nonzero_zero(ans, certanlyNonzeroIndex, maybeZeroIndex):
    """
    Find the indices of certainly non-zero and maybe zero values in the solution vector.

    Args:
        ans (numpy.ndarray): The solution vector.
        certanlyNonzeroIndex (list): The list of indices that are certainly non-zero.
        maybeZeroIndex (list): The list of indices that may be zero.
    """
    rest = [ans[i] for i in range(0, len(ans)) if i != certanlyNonzeroIndex[0] and i != maybeZeroIndex[0]]
    for i in range(0, len(rest)):
        tmp = np.max(np.abs(rest))
        tmp2 = int(np.argmax(np.abs(rest)))
        tmpind = find_val(rest[tmp2], ans)
        
        maxpos = sum([abs(ans[nonzind]) for nonzind in certanlyNonzeroIndex])/len(certanlyNonzeroIndex) - tmp
        minpos = tmp - abs(ans[maybeZeroIndex[0]])
        if maxpos < minpos:
            certanlyNonzeroIndex.append(tmpind)
        else:
            maybeZeroIndex.append(tmpind)
        
        rest = [ans[i] for i in range(0, len(ans)) if i not in certanlyNonzeroIndex and i not in maybeZeroIndex and i != tmpind]

#input parameters

#polynomyal stucture
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

#m
m_param = 4
#Нижній ліміт m
lowerbound_m = 1
#Верхній ліміт m
upperbound_m = 40
#n
input_size = 100

#for input vector
#Розмір вектора
inputvector_size = len(operations)
#Нижній ліміт вектора
inputvector_lowerbound = -10
#Верхній ліміт вектора
inputvector_upperbound = 10
#Кількість нулів в векторі
inutvector_zeros = 3

#Кількість повторень на експеримент
repeats_per_exp = 10
#Кількість експериментів
exp_count = 1000

#dispersion
dist_top = 10
d = [x for x in range(1, dist_top + 1)]

#Кількість експериментів для вибірки
exp_split = (repeats_per_exp // 10) * 6

#end of input parameters

err = np.zeros((dist_top, exp_count))
zeroes_wrong = np.zeros((dist_top, exp_count))

for f in range(dist_top):
    for k in range(exp_count):
        mat = generator.generate_input_matrix(input_size, m_param, lowerbound_m, upperbound_m)
        
        indices = generator.generate_input_vector(inputvector_size, inputvector_lowerbound, inputvector_upperbound, inutvector_zeros)
        y1, y2 = generator.generate_input_set(input_size, repeats_per_exp, d[f], indices, mat, exp_split, operations)
        yA1 = np.array([sum(y1[i])/exp_split for i in range(input_size)])
        
        ml1 = np.array([[operations[j](mat[i]) for j in range(inputvector_size)] for i in range(input_size)])
        
        ans = least_sqares(ml1, yA1)
        certanlyNonzeroIndex = [np.argmax(np.abs(ans))]
        maybeZeroIndex = [np.argmin(np.abs(ans))]
        
        find_nonzero_zero(ans, certanlyNonzeroIndex, maybeZeroIndex)
       
        comb = []
        for i in range(0, len(maybeZeroIndex)+ 1):
            comb.extend(list(combinations(maybeZeroIndex, i)))
            
        coef = calcBestInd(mat, certanlyNonzeroIndex, comb, y2,operations,yA1)
        ml2 = np.array([[operations[j](mat[i]) for j in range(inputvector_size) if j in coef] for i in range(input_size)])
        
        y = np.concatenate((y1, y2), axis=1)
        yA2 = np.array([sum(y[i])/(len(y[i])) for i in range(input_size)])
        ans2 = least_sqares(ml2, yA2)
        retind = []
        ansind = 0
        for i in range(len(indices)):
            if i in coef:
                retind.append(ans2[ansind])  
                ansind += 1
            else:
                retind.append(0)
        err[f, k],zeroes_wrong[f,k] = comMeasure(indices, retind)
        file = open("out.txt", "a")     
        
    data = np.zeros(5)
    data[0] = np.mean(err[f])
    data[1] = np.max(err[f])
    data[2] = np.min(err[f])
    data[3] = np.mean(zeroes_wrong[f])
    data[4] = np.max(zeroes_wrong[f])
    deploy_to_file("out.txt", data, d[f])
