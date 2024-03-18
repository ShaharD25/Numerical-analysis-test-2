#Date: 18.3.24
#Group members: Andrey Romanchuk ID;Shahar Dahan ID;Maya Rozenberg ID
#Source Git: https://github.com/lihiSabag/Numerical-Analysis-2023.git
#Private Git: https://github.com/ShaharD25/Numerical-analysis-test-2.git
#Name: shahar dahan

import numpy as np

def row_addition_elementary_matrix(n, target_row, source_row, scalar=1.0):
    if target_row < 0 or source_row < 0 or target_row >= n or source_row >= n:
        raise ValueError("Invalid row indices.")

    if target_row == source_row:
        raise ValueError("Source and target rows cannot be the same.")

    elementary_matrix = np.identity(n)
    elementary_matrix[target_row, source_row] = scalar

    return np.array(elementary_matrix)

def scalar_multiplication_elementary_matrix(n, row_index, scalar):

    if row_index < 0 or row_index >= n:
        raise ValueError("Invalid row index.")

    if scalar == 0:
        raise ValueError("Scalar cannot be zero for row multiplication.")

    elementary_matrix = np.identity(n)
    elementary_matrix[row_index, row_index] = scalar

    return np.array(elementary_matrix)

def solve_linear_system(coeff_matrix, b_vector):
    try:
        x_vector = np.linalg.solve(coeff_matrix, b_vector)
        return x_vector
    except np.linalg.LinAlgError:
        raise ValueError("System of equations is singular or not square, cannot be solved.")

def verify_inverse(original_matrix, inverse_matrix):
    product = np.dot(original_matrix, inverse_matrix)
    identity_matrix = np.identity(original_matrix.shape[0])

    # Check if the product is close to the identity matrix
    if np.allclose(product, identity_matrix):
        return True
    else:
        return False
def inverse(matrix):
    print(f"=================== Finding the inverse of a non-singular matrix using elementary row operations ===================\n {matrix}\n")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Input matrix must be square.")

    n = matrix.shape[0]
    identity = np.identity(n)

    # Perform row operations to transform the input matrix into the identity matrix
    for i in range(n):
        if matrix[i, i] == 0 and n == 2:
            raise ValueError("Matrix is singular, cannot find its inverse.")

        if matrix[i, i] == 0:
            matrix[i], matrix[i+1] = matrix[i+1].copy(), matrix[i].copy()
            print(f"The matrix after swap operation :\n {matrix}")
            print("------------------------------------------------------------------------------------------------------------------")


        if matrix[i, i] != 1:
            # Scale the current row to make the diagonal element 1
            scalar = 1.0 / matrix[i, i]
            elementary_matrix = scalar_multiplication_elementary_matrix(n, i, scalar)
            print(f"elementary matrix to make the diagonal element 1 :\n {elementary_matrix} \n")
            matrix = np.dot(elementary_matrix, matrix)
            print(f"The matrix after elementary operation :\n {matrix}")
            print("------------------------------------------------------------------------------------------------------------------")
            identity = np.dot(elementary_matrix, identity)

        # Zero out the elements above and below the diagonal
        for j in range(n):
            if i != j:
                scalar = -matrix[j, i]
                elementary_matrix = row_addition_elementary_matrix(n, j, i, scalar)
                print(f"elementary matrix for R{j+1} = R{j+1} + ({scalar}R{i+1}):\n {elementary_matrix} \n")
                matrix = np.dot(elementary_matrix, matrix)
                print(f"The matrix after elementary operation :\n {matrix}")
                print("------------------------------------------------------------------------------------------------------------------")
                identity = np.dot(elementary_matrix, identity)

    return identity



if __name__ == '__main__':

    A = np.array([[-1, 2, -3],
           [4, -5, 6],
           [7, 8, -9]])

    b = np.array([1, 2, 3])

    try:
        solution = solve_linear_system(A, b)
        A_inverse = inverse(A)
        print("\nSolution to the linear system: \n", solution)
        #print("\nInverse of matrix A: \n", A_inverse)
        #print("=====================================================================================================================")
        #if verify_inverse(A, A_inverse):
            #print("\nVerification: A * A_inverse is identity matrix.")
        #else:
            #print("\nVerification failed: A * A_inverse is not identity matrix.")

        #print("=====================================================================================================================")

    except ValueError as e:
        print(str(e))



