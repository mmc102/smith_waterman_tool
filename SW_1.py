import itertools
import numpy as np

string_a = 'GGTTGACTAGG'
string_b = 'GGTTTTTTTGACTA'


#form a matrix for string comparison 
def matrix(a,b,match_value= 3, gap_cost = 2):
    rix = np.zeros((len(a)+1,len(b)+1), np.int)
    for i,j in itertools.product(range(1,rix.shape[0]),range(1,rix.shape[1])):
        match = rix[i-1,j-1]+(match_value if a[i-1] == b[j-1] else - match_value)
        deletion = rix[i-1,j]-gap_cost
        insertion = rix[i,j - 1] - gap_cost
        rix[i,j] = max(match, deletion, insertion,0)
    return rix

def traceback(rix, b, b_='', old_i = 0):
    rix_flip = np.flip(np.flip(rix,0),1)
    i_,j_ = np.unravel_index(rix_flip.argmax(),rix_flip.shape)
    i,j = np.subtract(rix.shape,(i_+1,j_+1))
    if rix[i,j] == 0:
        return b_,j
    b_ = b[j-1]+ '-' +b_ if old_i - i > 1 else b[j-1] + b_
    return traceback(rix[0:i, 0:j], b, b_, i)

def smith_waterman(a,b,match_value = 3, gap_cost =2):
    a,b = a.upper(),b.upper()
    rix = matrix(a,b,match_value,gap_cost)
    b_,pos = traceback(rix,b)
    return pos,pos + len(b_)

print(matrix(string_a,string_b))

H = matrix(string_a,string_b)
print(traceback(H,string_b))

start,end = smith_waterman(string_a,string_b)

print(start,end)
print(string_a[start:end])