from hamiltonian import *

class ReducedDensity:
    
    def index(self,m1,m2):
        # m1 > m2
        return int(m1+self.S)*int(m1+self.S-1)/2+int(m2+self.S)
        
    def __init__(self,groundstate,hilbert):
        self.S = hilbert.S
        self.rank = (int(2*self.S)+1)*(int(2*self.S))/2
        self.N = hilbert.N
        self.reduced_density = numpy.zeros((self.rank,self.rank),dtype=complex)
        for index in range(hilbert.dim):
            if groundstate[index]!=0:
                vec = hilbert.search(index)
                for i in range(self.N):
                    for j in range(i+1,self.N):
                        m4,m3 = vec[i],vec[j] # m3>m4
                        self.reduced_density[self.index(m3,m4),self.index(m3,m4)]+=numpy.absolute(groundstate[index])**2
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
                                self.reduced_density[self.index(m2,m1),self.index(m3,m4)]+=sign*groundstate[index]*numpy.conjugate(groundstate[index2])
                                
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
                                self.reduced_density[self.index(m2,m1),self.index(m3,m4)]+=sign*groundstate[index]*numpy.conjugate(groundstate[index2])
                                
    def diagonalize(self):
        self.eigenvalues, self.eigenvectors = numpy.linalg.eigh(self.reduced_density)


            
if __name__=='__main__':
    
    S, N = 7.5, 6
    hilbert = HilbertSpace(S,N)
    
    pseudo_potential = [1]
    twobody = TwoBodyInteraction(pseudo_potential,S)
    onebody = OneBodyPotential([[0,0],[0.5,0.5,0.5]],S)
    
    hamiltonian = Hamiltonian(twobody, onebody, hilbert)
    hamiltonian.diagonalize(20)
    
    reduced_matrix = ReducedDensity(hamiltonian.groundstate, hilbert)
    reduced_matrix.diagonalize()
    print reduced_matrix.eigenvalues
    plt.figure()
    plt.plot([0]*reduced_matrix.rank, reduced_matrix.eigenvalues,marker = '_',ms=40)
    plt.ylim(-0.01,0.25)
    plt.show()
