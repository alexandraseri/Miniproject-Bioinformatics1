import sys


def isAMatch(x, y):
    if (x == 'A' and y == 'U') or (x == 'U' and y == 'A') or (x == 'G' and y == 'C') or (x == 'C' and y == 'G'):
        return 1

    else:
        return 0


if len(sys.argv) == 1:
    print("Error: No RNA sequence found.")

else:
    rna_string = sys.argv[1]
    length = len(rna_string)
    w, h = length, length
    matrix = [[0 for x in range(w)] for y in range(h)]
    for x in range(length-1):
        matrix[x][x+1] = isAMatch(rna_string[x], rna_string[x+1])

    for i in range(length-1):
        for j in range(i,length):
            matrix[i][j] = max(matrix[i+1][j-1]+isAMatch(rna_string[i], rna_string[j]), matrix[i+1][j], matrix[i][j-1])
            for k in range(i, j):
                matrix[i][j] = max(matrix[i][j], matrix[i][k] + matrix[k+1][j])

    for x in range(length - 1):
        print (matrix[x])