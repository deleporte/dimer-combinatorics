import scipy.optimize as optimization
import numpy as np
import scipy.linalg


def quad(x,a,b,c,d,e,f,g,h,i,j):
    val=a*x[0]*x[0]+b*x[0]*x[1]+c*x[1]*x[1]+d*x[2]*x[1]+e*x[0]*x[2]+f*x[2]*x[2]+g*x[0]+h*x[1]+i*x[2]+j
    return np.where(val>0,val,0)

def analyze(N=30):
    xdata=np.zeros((3,N**3)) #we will need at most N**3 points
    ydata=[]
    count=0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if exact_circle(N,i,j,k)>=1: #only consider the nonzero values
                    ydata.append(math.log(exact_circle(N,i,j,k))) #long int
                    xdata[0,count]=i
                    xdata[1,count]=j
                    xdata[2,count]=k
                    count += 1

    xdata=xdata.T[:count].T #trim the matrix
    p0=np.zeros(10)
    A,sigma = optimization.curve_fit(quad, xdata, np.array(ydata),p0)

    #eigenvalues and vectors of the 2nd order term
    mat=np.array([[A[0]*2,A[1],A[4]],[A[1],A[2]*2,A[3]],[A[4],A[3],2*A[5]]])
    values,vectors=scipy.linalg.eigh(mat)

    #maximal point
    vmax=np.dot(scipy.linalg.inv(mat),np.array([-A[6],-A[7],-A[8]]))
    dimmax,ptsmax,odimmax=vmax.tolist()

    #maximal value
    maxval=A[0]*dimmax*dimmax+A[1]*dimmax*ptsmax+A[2]*ptsmax*ptsmax+A[3]*ptsmax*odimmax+A[4]*dimmax*odimmax+A[5]*odimmax*odimmax+A[6]*dimmax+A[7]*ptsmax+A[8]*odimmax+A[9]
    
    return values,vectors,dimmax,ptsmax,odimmax,maxval



