import itertools
import numpy as np
import re

#strings, aligns string B to string A
A = 'GGTTGACTATT' 
B = 'TGTTCACTA'   

#substition matrix, this could be edited to penalize highly aginst unlikely subs
sub_array = np.array([[3,-3,-3,-3],[-3,3,-3,-3],[-3,-3,3,-3],[-3,-3,-3,3]])  #ATGC x ATGC
key = {'A':0,'T':1,'G':2,'C':3}
gap_cost = 2

#form the matrix
def make_matrix(a,b,sub_array,gap_cost):
    #init empty matrix
    M = np.zeros((len(a)+1,len(b)+1),np.int)
    #populate the matrix 
    for i,j in itertools.product(range(1,M.shape[0]),range(1,M.shape[1])):
        #generate match, insert and del value, take max and plop into matrix
        A_base = a[i-1]
        B_base = b[j-1]
        sub_cords = [(key.get(A_base)),(key.get(B_base))]
        match_value = sub_array[sub_cords[0]][sub_cords[1]]
        match = M[i-1,j-1]+(match_value if a[i-1] == b[j-1] else + match_value)  #if its a match, sum the score, else subtract 
        deletion = M[i-1,j]-gap_cost 
        insertion = M[i,j - 1]-gap_cost
        M[i,j] = max(match,deletion,insertion,0)
    print(M)
    return M

def traceback_matrix(M,b,b_ = '',old_i=0, cycle = 1, prev_i =0,prev_j=0, matches = 1,inserts = 0,deletions = 0,subs = 0):
    matrix_flip = np.flip(np.flip(M,0),1)
    i_,j_ = np.unravel_index(matrix_flip.argmax(),matrix_flip.shape)  #skips in a sub situation  
    i,j = np.subtract(M.shape,(i_+1,j_+1))

    #use i,j to trace the line, if we go left up or diagonal, can determine indels 
    if cycle == 1:
        matches += 1  #not sure why this needs to be here..added it because it makes the output correct...
        
    #what happens in a sub? not a match, but still goes diagonal, just regains some value after the sub
    if cycle > 1:
        if i+1 == prev_i and j+1 == prev_j:
            matches += 1
        elif i+2 == prev_i and j+2 == prev_j:
            print('*****SUB',b[j],'->',A[j])
            subs += 1
           
        elif i+1 == prev_i and j != prev_j:
            print('*****INS',b[j])
            inserts +=1
            b_ = ('('+b[j]+')') + b_

        elif i != prev_i and j+1 == prev_j:   
            print('*****DEL', A[j])
            deletions +=1
            
    #if bottom right of matrix is 0, return b_ = no string and j
    if M[i,j] == 0:
        return b_,j
    #otherwise, b_ contains the manipluations to b that best align with a
   
    b_ = b[j-1]+ '-' + b_ if old_i - i > 1 else b[j-1] + b_   #else handles insertions and matches, need conditon to prune the insertion for earmarking 
    print('mat:',matches)
    print('ins',inserts)
    print('dels',deletions)
    print('subs',subs)
    print('accuracy', (matches/(matches+inserts+deletions+subs))*100)
    print('------')
    cycle +=1
    #need to determine if we go diagonal (match), left (insertion) or right (deletion)
    
    return traceback_matrix(M[0:i, 0:j], b, b_, i,cycle, i, j,matches,inserts,deletions,subs)  #resizes the matrix, extends the aligned b,  old i becomes the latest i
   
def smith_waterman(a,b,sub_array, gap_cost):
    a,b = a.upper(),b.upper()
    M = make_matrix(a,b,sub_array,gap_cost)
    b_,pos  = traceback_matrix(M,b)
    b_stripped = re.sub("[\(\[].*?[\)\]]", "", b_)
    return pos, pos+(len(b_stripped)), b_




start, end, aligned = smith_waterman(A,B,sub_array,gap_cost)
print(A[start:end])
print(aligned)
#dont currently have differenation between subs and deletions, and it just doesnt give any indication of a sub. just doesnt incease the match