from hilbert import *
from potential import *
import numpy
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as plt 

               
class Hamiltonian:
        
    def __init__(self,twobody,onebody,hilbert):
        self.N = hilbert.N
        self.S = twobody.S
        self.dim = hilbert.dim
        row_indices = []
        index_pointer = [0]
        data = []
            
        for block in hilbert.vec_list:
            for vec in block:
                
                #diagonal elements need to be treated separately for csc form
                diagonal = 0
                for i in range(self.N):
                    for j in range(i+1,self.N):
                        m4,m3 = vec[i],vec[j]
                        diagonal+=twobody.matrix[int(m4+self.S),int(m3+self.S),int(m3+self.S),int(m4+self.S)]
                    diagonal+=onebody.matrix[int(vec[i]+self.S),int(vec[i]+self.S)]
                if (diagonal!=0):
                    row_indices.append(hilbert.search(vec))
                    data.append(diagonal)   
                
                #two-body interaction
                for i in range(self.N):
                    for j in range(i+1,self.N):
                        m4,m3 = vec[i],vec[j] # m3>m4                      
                        vec_proxy=list(vec)
                        sign = 1                    
                        m1,m2=m4,m3 # m2>m1
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
                                element = twobody.matrix[int(m1+self.S),int(m2+self.S),int(m3+self.S),int(m4+self.S)]
                                if(element!=0):
                                    row_indices.append(hilbert.search(tuple(sorted(vec_proxy))))
                                    data.append(sign*element)
                                
                        m1,m2=m4,m3
                        sign=1
                        while(m1>-self.S and m2<self.S):
                            m2+=1
                            m1-=1
                            if (m1 in vec or m2 in vec):
                                if (m1 in vec and m2 in vec):
                                    pass
                                else:
                                    sign*=-1                           
                            else:
                                vec_proxy[i],vec_proxy[j] = m1,m2
                                element = twobody.matrix[int(m1+self.S),int(m2+self.S),int(m3+self.S),int(m4+self.S)]
                                if(element!=0):
                                    row_indices.append(hilbert.search(tuple(sorted(vec_proxy))))
                                    data.append(sign*element)
                                
                #one-body
                for i in range(self.N):
                    m2 = vec[i]
                    vec_proxy = list(vec)
                    m1 = m2
                    sign = 1
                    while(m1>-self.S):
                        m1-=1
                        if m1 in vec:
                            sign*= -1
                        else:
                            vec_proxy[i]=m1
                            element = onebody.matrix[int(m1+self.S),int(m2+self.S)]
                            if (element!=0):
                                row_indices.append(hilbert.search(tuple(sorted(vec_proxy))))
                                data.append(sign*element)
                    
                    m1 = m2
                    sign = 1
                    while(m1<self.S):
                        m1+=1
                        if m1 in vec:
                            sign*= -1
                        else:
                            vec_proxy[i]=m1
                            element = onebody.matrix[int(m1+self.S),int(m2+self.S)]
                            if (element!=0):
                                row_indices.append(hilbert.search(tuple(sorted(vec_proxy))))
                                data.append(sign*onebody.matrix[int(m1+self.S),int(m2+self.S)])
                            
                index_pointer.append(len(data))
        self.matrix = csc_matrix((data,row_indices,index_pointer),shape=(self.dim,self.dim),dtype=complex)
                        
    def diagonalize(self,k):
        self.neig = k
        self.eigenvalues, self.eigenvectors = eigsh(self.matrix,k=k,which='SA')
        self.eigenvalues.sort()
        self.groundstate = self.eigenvectors[:,0]
        self.groundenergy = self.eigenvalues[0]   
        
    def print_groundstate(self,hilbert):
        for index in range(self.dim):
            if (numpy.absolute(self.groundstate[index])>1e-6):
                print hilbert.search(index), '\t', self.groundstate[index]
         
    def spectrum(self):
        self.distinct_eigenvalues = [self.eigenvalues[0]]
        angular_momentum = []
        prev = self.eigenvalues[0]
        prev_i=0
        for i in range(1,self.neig):
            if abs(self.eigenvalues[i]-prev)>1e-6:
                self.distinct_eigenvalues.append(self.eigenvalues[i])
                angular_momentum+=[(i-prev_i-1)/2]*(i-prev_i)
                prev = self.eigenvalues[i]
                prev_i = i
        return angular_momentum, self.eigenvalues[:prev_i]
        

        
if __name__ == '__main__':
    S, N = 4.5, 4
    hilbert = HilbertSpace(S,N)
    
    pseudo_potential = [1]
    twobody = TwoBodyInteraction(pseudo_potential,S)
    onebody = OneBodyPotential([],S)
    
    hamiltonian = Hamiltonian(twobody, onebody, hilbert)
    hamiltonian.diagonalize(100)
    #print hamiltonian.groundenergy
    angular_momentum, eigenvalues = hamiltonian.spectrum()
    print hamiltonian.eigenvalues
    print hamiltonian.distinct_eigenvalues
    '''
    plt.figure()
    plt.scatter(angular_momentum,eigenvalues)
    plt.show()
    '''