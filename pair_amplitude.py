from reducedmatrix import *
from randomstate import *

S,N = 4.5,4
hilbert = HilbertSpace(S,N)
zerohilbert = ZeroHilbert(hilbert)
zero_ang = ZeroAngularMomentum(zerohilbert)
zero_ang.diagonalize()
random_state = RandomState(zero_ang,10)

for state in random_state.state_list:
    reduced_matrix = ReducedDensity(state,zerohilbert)
    reduced_matrix.diagonalize()
    eigenvalues = reduced_matrix.eigenvalues
    distinct_eigv = [0]*len(range(int(2*S)-1,-1,-2))
    prev_i = 0
    prev = eigenvalues[0]
    for i in range(1,reduced_matrix.rank):
        if (abs(eigenvalues[i]-eigenvalues[prev_i])>1e-14):
            l = (i-prev_i-1)/2
            distinct_eigv[(int(2*S)-1-l)/2]=prev
            prev_i = i
            prev = eigenvalues[i]
    i+=1
    distinct_eigv[(int(2*S)-1-(i-prev_i-1)/2)/2]=prev
    print distinct_eigv