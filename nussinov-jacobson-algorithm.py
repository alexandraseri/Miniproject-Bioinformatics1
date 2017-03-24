import re
import sys



# Checks if two bases can make a bp by Watson-Crick
# @param: x,y the bases to check
# @return: 1 if are a bp, 0 otherwise
def isAMatch(x, y):
    if (x == 'A' and y == 'U') or (x == 'U' and y == 'A') or (x == 'G' and y == 'C') or (x == 'C' and y == 'G'):
        return 1

    else:
        return 0


# Checks program args for problems with input (No proper RNA sequence)
# return: True if there is a proper RNA sequence, False otherwise
def checkArgs():
    if len(sys.argv) < 2: # no string input found
        print("Error: No RNA sequence found.")
        return False
    
    seq = sys.argv[1]
    if bool(re.match('^[AUCG]*$',seq)) is False: # the string is not a RNA sequence
        print("Error: Not a RNA sequence.")
        return False
    
    return True


# Calculates the matrix for maximum number of bp for secondery structure.abs
# @param: seq the RNA sequence to analyze.
# @return: the matrix after it's calculated.
def calculateMatrix(seq):
    length = len(seq)
    w, h = length, length
    matrix = [[0 for x in range(w)] for y in range(h)]  # matrix initialization

    for j in range(2,length): # the distance between 2 bases
        for i in range(length-j): # the bases are i and i+j
            match = isAMatch(seq[i], seq[i+j])
            if (match == 1): # if there is a match between bases, save the curren number of matches so far in 'match' var
                match = matrix[i+1][i+j-1]+match
                
            temp = 0
            for k in range(i+1, i+j): # calculate for i<k<j if there is a better splicing in sequence for matching
                temp = max(temp, matrix[i][k] + matrix[k+1][i+j])
                
            matrix[i][i+j] = max(match,temp, matrix[i+1][i+j], matrix[i][i+j-1]) # max between 4 options from dp formula
    
    return matrix


# Traceback for 'dot bracket notation' from matrix
# @param: matrix the matrix to traceback from
# @param: seq the RNA sequence to traceback from
# @return: the answer in a 'dot bracket notation' form.
def traceback(matrix,seq):
    length = len(seq)
    answer = ['.' for x in range(length)] # answer initialization with no bps

    stack = []
    stack.append((0,length-1))

    while(len(stack)!= 0): # until there are no more bps
        pair = stack.pop()
        i,j = pair[0],pair[1] 
        if(i<j): # check only if i<j
            if (isAMatch(seq[i], seq[j]) == 1 and matrix[i+1][j-1] + 1 == matrix[i][j]): # match between i, j
                answer[i] = '('
                answer[j] = ')'
                stack.append((i+1,j-1))
                
            elif (matrix[i][j-1] == matrix[i][j]): # match between i, j-1
                stack.append((i,j-1))
                
            elif(matrix[i+1][j] == matrix[i][j]): # match between i+1,j
                stack.append((i+1,j))
                    
            else:
                for k in range(i+1,j):
                    if(matrix[i][k]+matrix[k+1][j] == matrix[i][j]): #match between i,k and k+1,j
                        stack.append((k+1,j))
                        stack.append((i,k))
                        break

    return ''.join(answer)


def main():
    if checkArgs() is False: # check program args
        return
        
    else:
        seq = sys.argv[1]
        matrix = calculateMatrix(seq) # calculate matrix
        answer = traceback(matrix, seq) # calculate answer
        
        print(answer)
        return

if __name__ == "__main__":
    main()