from zerohilbert import *
import numpy
import math
import sys

class ZeroAngularMomentum:
    
    def __init__(self, zerohilbert):
        self.S = zerohilbert.S
        self.N = zerohilbert.N
        self.rank = zerohilbert.dim
        self.matrix = numpy.zeros((self.rank,self.rank))
        for vec in zerohilbert.vec_list:
            index1 = zerohilbert.search(vec)
            self.matrix[index1,index1]+=self.N*self.S*(self.S+1)  
            for i in range(self.N):
                for j in range(i+1,self.N):
                    vec_proxy = list(vec)
                    self.matrix[index1,index1]+=2*vec[i]*vec[j]
                    if (vec[i]+1==vec[j]):
                        self.matrix[index1,index1]-=self.S*(self.S+1)-vec[i]*vec[j]
                    if (vec[j]>vec[i]+2 and vec[i]+1 not in vec and vec[j]-1 not in vec):
                        vec_proxy[i] = vec[i]+1
                        vec_proxy[j] = vec[j]-1
                        self.matrix[zerohilbert.search(tuple(vec_proxy)),index1]+=\
                        math.sqrt(self.S*(self.S+1)-vec[i]*(vec[i]+1))*math.sqrt(self.S*(self.S+1)-vec[j]*(vec[j]-1))
                    if (vec[i]>-self.S and vec[j]<self.S and vec[i]-1 not in vec and vec[j]+1 not in vec):
                        vec_proxy[i] = vec[i]-1
                        vec_proxy[j] = vec[j]+1
                        self.matrix[zerohilbert.search(tuple(vec_proxy)),index1]+=\
                        math.sqrt(self.S*(self.S+1)-vec[i]*(vec[i]-1))*math.sqrt(self.S*(self.S+1)-vec[j]*(vec[j]+1))
                        
    def diagonalize(self):
        self.eigenvalues, self.eigenvectors = numpy.linalg.eigh(self.matrix)
        self.zero_eigvec = []
        self.nzero = 0
        while abs(self.eigenvalues[self.nzero])<1e-10: self.nzero+=1
        self.zero_eigvec = self.eigenvectors[:,:self.nzero]
                
                        
if __name__ == '__main__':
    S,N = 4.5,4
    hilbert = HilbertSpace(S,N)
    zerohilbert = ZeroHilbert(hilbert)
    zero_ang = ZeroAngularMomentum(zerohilbert)
    zero_ang.diagonalize()
    print zero_ang.nzero
    print zero_ang.zero_eigvec