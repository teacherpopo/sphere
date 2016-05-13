import random
import math
from zeroam import *

class RandomState:
    
    def random_generator(self,lentuple,ntuple,real):
        self.random_list = []
        for i in range(ntuple):
            group = ()
            square_sum = 0
            for j in range(lentuple):
                if real==False:
                    group+=(complex(random.random(),random.random()),)
                else:
                    group+=(complex(random.random(),0),)
                square_sum+=numpy.absolute(group[-1])**2
            group = map(lambda x:x/math.sqrt(square_sum),group)
            self.random_list.append(group)
        
        
    def __init__(self,zeroam,nstate,real=False):
        self.nzero = zeroam.nzero
        self.rank = zeroam.rank
        self.random_generator(self.nzero,nstate,real)
        self.state_list=[numpy.dot(zeroam.zero_eigvec,group) for group in self.random_list]
                
if __name__ == '__main__':
    S,N = 4.5,4
    hilbert = HilbertSpace(S,N)
    zerohilbert = ZeroHilbert(hilbert)
    zero_ang = ZeroAngularMomentum(zerohilbert)
    zero_ang.diagonalize()
    random_state = RandomState(zero_ang,10)
    print random_state.random_list
    print random_state.state_list
