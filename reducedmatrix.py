from __future__ import print_function
from hamiltonian import *
import numpy


class ReducedDensity:
    
    def index(self,m1,m2):
        # m1 > m2
        return int(m1+self.S)*int(m1+self.S-1)/2+int(m2+self.S)
        
    def lm_index(self,l,m):
        index = 0
        for i in range(int(2*self.S)-1,l,-2):
            index+=2*i+1
        return index+m+l
        
    def search(self,index):
        block_index = 0
        i = int(2*self.S)-1
        while(block_index+2*i+1<=index):
            block_index+=2*i+1
            i-=2
        return (i,index-block_index-i)
        
    def __init__(self,groundstate,hilbert):
        self.S = hilbert.S
        self.rank = (int(2*self.S)+1)*(int(2*self.S))/2
        self.N = hilbert.N
        self.matrix = numpy.zeros((self.rank,self.rank),dtype=complex)
        for index in range(hilbert.dim):
            if groundstate[index]!=0:
                vec = hilbert.search(index)
                for i in range(self.N):
                    for j in range(i+1,self.N):
                        m4,m3 = vec[i],vec[j] # m3>m4
                        self.matrix[self.index(m3,m4),self.index(m3,m4)]+=numpy.absolute(groundstate[index])**2
                        vec_proxy=list(vec)
                        sign = 1                    
                        m1,m2=m4,m3
                        while(m2>m1+2):
                            m2-=1
                            m1+=1
                            if (m1 in vec or m2 in vec):
                                if (m1 in vec and m2 in vec):
                                    pass
                                else:
                                    sign*=-1                           
                            else:
                                vec_proxy[i],vec_proxy[j] = m1,m2
                                index2 = hilbert.search(tuple(sorted(vec_proxy)))
                                self.matrix[self.index(m2,m1),self.index(m3,m4)]+=sign*groundstate[index]*numpy.conjugate(groundstate[index2])
                                
                        m1,m2=m4,m3
                        sign=1
                        while(m1>-S and m2<S):
                            m2+=1
                            m1-=1
                            if (m1 in vec or m2 in vec):
                                if (m1 in vec and m2 in vec):
                                    pass
                                else:
                                    sign*=-1                           
                            else:
                                vec_proxy[i],vec_proxy[j] = m1,m2
                                index2 = hilbert.search(tuple(sorted(vec_proxy)))
                                self.matrix[self.index(m2,m1),self.index(m3,m4)]+=sign*groundstate[index]*numpy.conjugate(groundstate[index2])
                                
    def basis_transform(self):
        self.basis = numpy.zeros((self.rank,self.rank))
        for l in range(int(2*self.S)-1,-1,-2):
            for m in range(-l,l+1):
                for m1 in numpy.arange(-self.S,self.S+1,1):
                    m2 = m-m1
                    if (m2>m1 and m2<=self.S):
                        self.basis[self.lm_index(l,m),self.index(m2,m1)]+=(numpy.sqrt(2)*sympy.N(CG(S,m2,S,m1,l,m).doit()))
                        
        self.lm_matrix = numpy.dot(numpy.dot(self.basis, self.matrix),numpy.transpose(self.basis))
                                
    def diagonalize(self):
        self.eigenvalues, self.eigenvectors = numpy.linalg.eigh(self.lm_matrix)
        
    def print_pair(self):
        for i in range(self.rank):
            if (i>0 and abs(self.eigenvalues[i]-self.eigenvalues[i-1])>1e-14): print()
            print('%f: '%self.eigenvalues[i],end='')
            for index in range(self.rank):
                modulus = numpy.absolute(self.eigenvectors[index,i])
                if (modulus>1e-14):
                    print('%f (%d,%d) '%((modulus,)+self.search(index)),end='')
            print()
            


            
if __name__=='__main__':
    
    S, N = 3, 3
    hilbert = HilbertSpace(S,N)
    
    pseudo_potential = [1]
    onebody_potential = [[0,0.1]]
    twobody = TwoBodyInteraction(pseudo_potential,S)
    onebody = OneBodyPotential(onebody_potential,S)
    
    hamiltonian = Hamiltonian(twobody, onebody, hilbert)
    hamiltonian.diagonalize(1)
    
    reduced_matrix = ReducedDensity(hamiltonian.groundstate, hilbert)
    reduced_matrix.basis_transform()
    reduced_matrix.diagonalize()
    reduced_matrix.print_pair()
    '''
    plt.figure()
    plt.plot([0]*reduced_matrix.rank, reduced_matrix.eigenvalues,marker = '_',ms=40)
    plt.ylim(-0.01,0.25)
    plt.show()'''
