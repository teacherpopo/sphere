from reducedmatrix import *
from randomstate import *
from scipy.spatial import ConvexHull
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

S,N = 4.5,4
hilbert = HilbertSpace(S,N)
zerohilbert = ZeroHilbert(hilbert)
zero_ang = ZeroAngularMomentum(zerohilbert)
zero_ang.diagonalize()
for iteration in range(1):
    random_state = RandomState(zero_ang,500,real=False)
    
    points = []
    for state in random_state.state_list:
        reduced_matrix = ReducedDensity(state,zerohilbert)
        reduced_matrix.diagonalize()
        eigenvalues = reduced_matrix.eigenvalues
        distinct_eigv = [0]*len(range(int(2*S)-1,-1,-2))
        prev_i = 0
        prev = eigenvalues[0]
        for i in range(1,reduced_matrix.rank):
            if (abs(eigenvalues[i]-eigenvalues[prev_i])>1e-10):
                l = (i-prev_i-1)/2
                distinct_eigv[(int(2*S)-1-l)/2]=prev
                prev_i = i
                prev = eigenvalues[i]
        i+=1
        l = (i-prev_i-1)/2
        distinct_eigv[(int(2*S)-1-l)/2]=prev
        points.append(distinct_eigv)
    #print zero_ang.nzero
    #print len(distinct_eigv)   
    #hull = ConvexHull(points)
    if(iteration==0):variance, principle = numpy.linalg.eigh(numpy.cov(numpy.transpose(points)))
    n = 0
    while (abs(variance[n])<1e-10):n+=1
    x,y,z,w,v = zip(*points)
    #fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    #ax.scatter(x,y,z)
    #print len(distinct_eigv)-n # 2 5 2, 2 6 2, 6 8 6, 10 9 7
    new_points = []
    for point in points:
        new_point=()
        for i in range(n,len(distinct_eigv)):
            new_point+=(numpy.dot(principle[:,i],point),)
        new_points.append(new_point)
    #plt.scatter(*zip(*new_points),color='r')
    #plt.show()
