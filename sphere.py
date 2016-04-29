from sympy.physics.quantum.cg import CG
from scipy.misc import comb
import numpy
import sympy
import sys
import matplotlib.pyplot as plt 

class HilbertSpace:
    
    def buildspace(self,vec_list,vec,S,N,start):
        if (N==1):
            for i in range(start,S+1):
                vec_list[sum(vec+(i,))-self.lmin].append(vec+(i,))
        else:
            for i in range(start,S+1):
                self.buildspace(vec_list,vec+(i,),S,N-1,i+1)
                
    def __init__(self,S,N):
        self.lmin = sum(range(-S,-S+N))
        self.lmax = -self.lmin
        self.vec_list = [[] for _ in range((self.lmax-self.lmin+1))]
        self.buildspace(self.vec_list,(),S,N,-S)
        
class TwoBodyInteraction:
    
    def __init__(self,pseudo_potential,S):
        self.S=S
        self.rank = 2*S+1
        self.twobody=numpy.zeros((self.rank,self.rank,self.rank,self.rank))
        for i in range(len(pseudo_potential)):
            self.l = 2*i+1
            self.potential = pseudo_potential[i]
            for m in range(-(2*S-self.l),2*S-self.l+1):
                for m1 in range(-S,S+1):
                    self.m2 = m-m1
                    if (self.m2>m1 and self.m2<=S):
                        for m3 in range(-S,S+1):
                            self.m4 = m-m3
                            if(self.m4<m3 and self.m4>=-S):
                                self.twobody[m1,self.m2,m3,self.m4]+=-sympy.N(CG(S,m1,S,self.m2,2*S-self.l,m).doit())\
                                *sympy.N(CG(S,m3,S,self.m4,2*S-self.l,m).doit())*self.potential*4
                                
    def __str__(self):
        self.str=''
        for m1 in range(-self.S,self.S+1):
            for m2 in range(m1+1,self.S+1):        
                for m4 in range(-self.S,self.S+1):
                    for m3 in range(m4+1,self.S+1):
                        if (self.twobody[m1,m2,m3,m4]!=0):
                            self.str += ','.join([str(x) for x in [m1,m2,m3,m4]])+': '+str(self.twobody[m1,m2,m3,m4])+'\n'
        return self.str

def search(vec,vec_list):
        tail = len(vec_list)-1
        head = 0
        i = (head+tail)/2
        while(vec!=vec_list[i]):
            if vec<vec_list[i]:
                tail = i-1
                i=(head+tail)/2
            else:
                head = i+1
                i=(head+tail)/2
        return i
               
class Hamiltonian:
        
    def __init__(self,twobody,hilbert,S,N):
        self.hamil_list = []
        for block in hilbert.vec_list:
            dim = len(block)
            hamiltonian = numpy.zeros((dim,dim))
            for index1 in range(dim):
                vec = block[index1]
                for i in range(N):
                    for j in range(i+1,N):
                            else:
                                sign*=-1
                        
                        else:
                            vec_proxy[i],vec_proxy[j] = m1,m2
                            index2 = search(tuple(sorted(vec_proxy)),hilbert.vec_list)
                            self.hamiltonian[index2,index1]+=sign*twobody.twobody[m1,m2,m3,m4]
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
                            index2 = search(tuple(sorted(vec_proxy)),hilbert.vec_list)
                            self.hamiltonian[index2,index1]+=sign*twobody.twobody[m1,m2,m3,m4]
     
    def print_matrix(self):
        for i in range(self.dim):
            for j in range(self.dim):
                sys.stdout.write(str(self.hamiltonian[i,j])+', ')
            print '\n'
                        
    def diagonalize(self):
        
    def print_groundstate(self,hilbert):
        for index in range(self.dim):
            if (abs(self.groundstate[index])>1e-6):
                print hilbert.vec_list[index], '\t', self.groundstate[index]
         
    def plot_spectrum(self):
        angular_momentum = []
        eigen_sorted = sorted(self.eigenvalues)
        prev = eigen_sorted[0]
        prev_i=0
        for i in range(1,len(eigen_sorted)):
            if abs(eigen_sorted[i]-prev)>1e-6:
                angular_momentum+=[(i-prev_i-1)/2]*(i-prev_i)
                prev = eigen_sorted[i]
                prev_i = i
        i+=1
        angular_momentum+=[(i-prev_i-1)/2]*(i-prev_i)
        plt.scatter(angular_momentum, eigen_sorted)

class ReducedDensity:
    def __init__(self,groundstate,hilbert,S):
        self.dim = 2*S+1
        self.N = len(hilbert.vec_list[0])
        self.reduced_density = numpy.zeros((self.dim**2,self.dim**2))
        for index in range(len(hilbert.vec_list)):
            if groundstate[index] != 0.0:
                vec = hilbert.vec_list[index]
                for i in range(N):
                    for j in range(i+1,N):
                        m4,m3 = vec[i],vec[j] # m3>m4
                        self.reduced_density[(m3+S)*self.dim+m4+S,(m3+S)*self.dim+m4+S]+=groundstate[index]**2
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
                                index2 = search(tuple(sorted(vec_proxy)),hilbert.vec_list)
                                self.reduced_density[(m2+S)*self.dim+m1+S,(m3+S)*self.dim+m4+S]+=sign*groundstate[index]*groundstate[index2]
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
                                index2 = search(tuple(sorted(vec_proxy)),hilbert.vec_list)
                                self.reduced_density[(m2+S)*self.dim+m1+S,(m3+S)*self.dim+m4+S]+=sign*groundstate[index]*groundstate[index2]
                                
    def diagonalize(self):
        self.eigenvalues, self.eigenvectors = numpy.linalg.eigh(self.reduced_density)
        
                
        
            
if __name__=='__main__':
    
    hilbertspace = HilbertSpace(S,N)
    
    twobody = TwoBodyInteraction([1,0.1,0.01],S)
    #print twobody
    
    hamiltonian = Hamiltonian(twobody, hilbertspace, S,N)
    hamiltonian.diagonalize()
    #print hamiltonian.eigenvalues
    #print hamiltonian.eigenvectors
    #hamiltonian.print_groundstate(hilbertspace)
    plt.show()