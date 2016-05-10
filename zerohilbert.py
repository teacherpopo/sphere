from hilbert import *

class ZeroHilbert:
    
    def __init__(self,hilbert):
        if (int(2*hilbert.lmax)%2==1):
            print 'Not compatible with Laughlin state'
            return
        self.vec_list = hilbert.vec_list[int(hilbert.lmax)]
        self.dim = len(self.vec_list)
        self.S = hilbert.S
        self.N = hilbert.N
        
    def search(self,vec):
        if (type(vec) == tuple):
            tail = self.dim-1
            head = 0
            i = (head+tail)/2
            while(vec!=self.vec_list[i]):
                if vec<self.vec_list[i]:
                    tail = i-1
                else:
                    head = i+1
                i=(head+tail)/2
            return i
        else:
            return self.vec_list[vec]
            
if __name__ == '__main__':
    S,N = 3,3
    hilbert = HilbertSpace(S,N)
    zerohilbert = ZeroHilbert(hilbert)
    print zerohilbert.vec_list