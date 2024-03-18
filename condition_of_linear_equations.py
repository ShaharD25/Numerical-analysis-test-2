#Date: 18.3.24
#Group members: Andrey Romanchuk ID;Shahar Dahan ID;Maya Rozenberg ID
#Source Git: https://github.com/lihiSabag/Numerical-Analysis-2023.git
#Private Git:https://github.com/ShaharD25/Numerical-analysis-test-2.git
#Name: shahar dahan

import numpy as np
from inverse_matrix import inverse

def print_matrix(matrix):
    for row in matrix:
        for element in row:
            print(element, end=" ")  # Print each element in the row
        print()  # Move to the next row
    print()



def norm(mat):
    size = len(mat)
    max_row = 0
    for row in range(size):
        sum_row = 0
        for col in range(size):
            sum_row += abs(mat[row][col])
        if sum_row > max_row:
            max_row = sum_row
    return max_row


def condition_number(A):
    # Step 1: Calculate the max norm (infinity norm) of A
    norm_A = norm(A)   #works great

    # Step 2: Calculate the inverse of A
    A_inv = inverse(A)

    # Step 3: Calculate the max norm of the inverse of A
    norm_A_inv = norm(A_inv)

    # Step 4: Compute the condition number
    cond = norm_A * norm_A_inv

    print("A:")
    print_matrix(A)

    print("inverse of A:")
    print_matrix(A_inv)

    print("Max Norm of A:", norm_A, "\n")

    print("max norm of the inverse of A:", norm_A_inv)

    return cond


if __name__ == '__main__':
    A = np.array([[2, 2, 3],
                  [4, -6, 2],
                  [4,  7, 5]])


    cond = condition_number(A)

    print("\n condition number: ", cond)





