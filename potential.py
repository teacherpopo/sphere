from sympy.physics.quantum.cg import CG
import numpy
import sympy

class TwoBodyInteraction:
    
    def __init__(self,pseudo_potential,S):
        self.S=S
        self.rank = int(2*S+1)
        self.matrix=numpy.zeros((self.rank,self.rank,self.rank,self.rank))
        for i in range(len(pseudo_potential)):
            l = 2*i+1
            potential = pseudo_potential[i]
            for m in range(-(int(2*S)-l),int(2*S)-l+1):
                for m1 in numpy.arange(-S,S+1,1):
                    m2 = m-m1
                    if (m2>m1 and m2<=S):
                        for m3 in numpy.arange(-S,S+1,1):
                            m4 = m-m3
                            if(m4<m3 and m4>=-S):
                                self.matrix[int(m1+S),int(m2+S),int(m3+S),int(m4+S)]+=\
                                -sympy.N(CG(S,m1,S,m2,2*S-l,m).doit())\
                                *sympy.N(CG(S,m3,S,m4,2*S-l,m).doit())*potential*4 # 4 because we only consider m2>m1 and m4<m3
                                
    def __str__(self):
        string='Two-body Interaction:\n'
        for m1 in numpy.arange(-self.S,self.S+1):
            for m2 in numpy.arange(m1+1,self.S+1):        
                for m4 in numpy.arange(-self.S,self.S+1):
                    for m3 in numpy.arange(m4+1,self.S+1):
                        temp = self.matrix[int(m1+self.S),int(m2+self.S),int(m3+self.S),int(m4+self.S)]
                        if (temp!=0):
                            string += ','.join([str(x) for x in [m1,m2,m3,m4]])+': '+str(temp)+'\n'
        return string
        
class OneBodyPotential:
    def __init__(self, potential_list, S):
        self.S = S
        self.rank = int(2*S)+1
        self.matrix = numpy.zeros((self.rank,self.rank),dtype=complex)
        for i in range(len(potential_list)):
            l=i+1
            potential_l = potential_list[i]
            for m in range(l+1):
                potential = potential_l[m]
                for m2 in numpy.arange(-S,S+1,1):
                    m1 = m+m2
                    if(abs(m1)<=S):
                        self.matrix[int(m1+S),int(m2+S)]+=potential\
                        *sympy.N(CG(S,m2,l,m,S,m1).doit())
                if (m!=0):
                    potential = (-1)**m*numpy.conjugate(potential)
                    m = -m
                    for m2 in numpy.arange(-S,S+1,1):
                        m1 = m+m2
                        if(abs(m1)<=S):
                            self.matrix[int(m1+S),int(m2+S)]+=potential\
                            *sympy.N(CG(S,m2,l,m,S,m1).doit())
     
    def __str__(self):
        string = 'One-body potential:\n'
        for m1 in numpy.arange(-self.S,self.S+1,1):
            for m2 in numpy.arange(-self.S,self.S+1,1):
                string+=str(self.matrix[int(m1+self.S),int(m2+self.S)])+'\t'
            string+='\n'
        return string      
            

if __name__ == '__main__':
    S = 1.5
    pseudo_potential = [0,0.1]
    twobody = TwoBodyInteraction(pseudo_potential,S)
    onebody = OneBodyPotential([[1,0.1]],S)
    print twobody
    print onebody