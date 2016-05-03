import numpy

class HilbertSpace:
    def search(self,vec):
        if (type(vec) == tuple):
            block_index = int(sum(vec)-self.lmin)
            block = self.vec_list[block_index]
            tail = len(block)-1
            head = 0
            i = (head+tail)/2
            while(vec!=block[i]):
                if vec<block[i]:
                    tail = i-1
                    i=(head+tail)/2
                else:
                    head = i+1
                    i=(head+tail)/2
            return self.block_start[block_index]+i
        else:
            i=0
            while(vec>=self.block_start[i+1]): i+=1
            return self.vec_list[i][vec-self.block_start[i]]
    
    def buildspace(self,vec_list,vec,S,N,start):
        if (N==1):
            for i in numpy.arange(start,S+1,1):
                vec_list[int(sum(vec+(i,))-self.lmin)].append(vec+(i,))
        else:
            for i in numpy.arange(start,S+1,1):
                self.buildspace(vec_list,vec+(i,),S,N-1,i+1)
                
    def __init__(self,S,N):
        self.S = S
        self.N = N
        self.lmin = sum([-S+x for x in range(N)])
        self.lmax = -self.lmin
        self.vec_list = [[] for _ in range(int(self.lmax-self.lmin+1))]
        self.buildspace(self.vec_list,(),S,N,-S)
        
        self.block_start = [0]       
        for block in self.vec_list:
            self.block_start.append(self.block_start[-1]+len(block))
        self.dim = self.block_start[-1]
        
if __name__ == '__main__':
    S, N = 2.5, 3
    hilbert = HilbertSpace(S,N)
    print hilbert.vec_list
    
    
    