#Date: 18.3.24
#Group members: Andrey Romanchuk ID;Shahar Dahan ID;Maya Rozenberg ID
#Source Git: https://github.com/lihiSabag/Numerical-Analysis-2023.git
#Private Git: https://github.com/ShaharD25/Numerical-analysis-test-2.git
#Name: shahar dahan

def isDDM(matrix, n):
    for i in range(0, n):
        sum = 0
        for j in range(0, n):
            sum = sum + abs(matrix[i][j])
        sum = sum - abs(matrix[i][i])
        if (abs(matrix[i][i]) < sum):
            return False

    return True

def Determinant(matrix, mul):
    width = len(matrix)
    if width == 1:
        return mul * matrix[0][0]
    else:
        sign = -1
        det = 0
        for i in range(width):
            m = []
            for j in range(1, width):
                buff = []
                for k in range(width):
                    if k != i:
                        buff.append(matrix[j][k])
                m.append(buff)
            sign *= -1
            det = det + mul * Determinant(m, sign * matrix[0][i])
    return det


def PivotMtx(matrix, vector):
    for i in range(len(matrix)):
        max = matrix[i][i]
        flag = i
        for j in range(i, len(matrix)):
            if(matrix[i][j] > max):
                max = matrix[i][j]
                flag = j
        if(flag != i):
            matrix[i], matrix[j] = matrix[j], matrix[i]
            vector[i], vector[j] = vector[j], vector[i]

    return matrix, vector


def MakeIMatrix(cols, rows):
    return [[1 if x == y else 0 for y in range(cols)] for x in range(rows)]


def findLDU(matrix):
    L, D, U = MakeIMatrix(len(matrix), len(matrix)), MakeIMatrix(len(matrix), len(matrix)), MakeIMatrix(len(matrix), len(matrix))

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i > j:
                L[i][j] = matrix[i][j]
            elif i == j:
                L[i][i], U[i][i], D[i][i] = 0, 0, matrix[i][i]
            else:
                U[i][j] = matrix[i][j]

    return L, D, U


def findGH(matrix, k):
    L, D, U = findLDU(matrix)
    H = getMatrixInverse(D)
    if k == "Jacobi":
        H = getMatrixInverse(D)
    else:
        H = getMatrixInverse(addmatrix(L, D))
    G = copy_matrix(H)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            G[i][j] = -G[i][j]
    if(k == "Jacobi"):
        G = matrix_multiply(G, addmatrix(L, U))
    else:
        G = matrix_multiply(G, U)

    return G, H

def addmatrix(matrixA, matrixB):
    return [[matrixA[i][j] + matrixB[i][j] for j in range
(len(matrixA[0]))] for i in range(len(matrixA))]


def copy_matrix(matrix):
    rows = len(matrix)
    cols = len(matrix[0])

    MC = [[0 for _ in range(len(matrix))] for _ in range(len(matrix))]

    for i in range(rows):
        for j in range(cols):
            MC[i][j] = matrix[i][j]

    return MC

def matrix_multiply(matrixA, matrixB):
    matrixC = [[0.0] * len(matrixB[0]) for _ in range(len(matrixA))]
    for i in range(len(matrixA)):
        for j in range(len(matrixB[0])):
            for k in range(len(matrixB)):
                matrixC[i][j] = matrixC[i][j] + matrixA[i][k] * matrixB[k][j]

    return matrixC


def transposeMatrix(matrix):
    return list(map(list, zip(*matrix)))

def getMatrixMinor(matrix, i, j):
    return [row[:j] + row[j+1:] for row in (matrix[:i] + matrix[i + 1:])]

def getMatrixInverse(matrix):
    determinant = Determinant(matrix, 1)
    if len(matrix) == 2:
        return [[matrix[1][1] / determinant, -1 * matrix[0][1] / determinant],
                [-1 * matrix[1][0] / determinant, matrix[0][0] / determinant]]

    cofactors = []
    for r in range(len(matrix)):
        cofactorRow = []
        for c in range(len(matrix)):
            minor = getMatrixMinor(matrix, r, c)
            cofactorRow.append(((-1)**(r+c)) * Determinant(minor, 1))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors


def JacobiandGausseSeidel(matrix, vector, method):
    matrix, vector = PivotMtx(matrix, vector)
    count = -1
    numberofint = 11

    if isDDM(matrix, len(matrix)) == False:
        print("Its not diagonally dominant Matrix")
        count = 150
    if method == "Jacobi":
        G, H = findGH(matrix, "Jacobi")
    else:
        G, H = findGH(matrix, "GausseSeidel")
    tmp = [[0 for _ in range(1)] for _ in range(len(vector))]
    print("0." + str(tmp))
    while count != 0:
        new_vector = tmp
        tmp = addmatrix(matrix_multiply(G, tmp), matrix_multiply(H, vector))
        print(str(numberofint) + "." + str(tmp) + "\n")
        flag = 0
        for i in range(len(tmp)):
            if abs(new_vector[i][0] - tmp[i][0]) >= 0.001:
                flag = 1
        if flag == 0:
            break
        count -= 1
        numberofint += 1
    return tmp


# Driver Code

matrixA = [[-1, 1, 3, -3, 1],
           [3, -3, -4, 2, 3],
           [2, 1, -5, -3, 5],
           [-5, -6, 4, 1, 3],
           [3, -2, -2, -3, 5]]
vectorB = [[3], [8], [2], [14], [6]]

while(True):
    print("Welcome to EX2 - Jaacobi/Gausse-Seidel")
    print("Please choose ")
    selection = input("[1] Jaacobi OR [2] Gausse-Seidel (By pressing other number the program will shutdown.): ")
    if selection == "1":
        JacobiandGausseSeidel(matrixA, vectorB, "Jacobi")
    elif selection == "2":
        JacobiandGausseSeidel(matrixA, vectorB, "GausseSeidel")
    else:
        break
