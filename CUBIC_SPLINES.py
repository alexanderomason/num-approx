
import numpy as np



def spln(X,F):
    '''where X is a list of nodes
    and F is a list of function values at those nodes,
    the output is a list two two entries [M,G],
    M is the list-list-matrix needed to solve for the z_js
    G is the list-vector used to create the system '''
    n=len(X)-1
    H=[]
    for k in range(n):
        H.append(X[k+1]-X[k])
    M=Zero(n-1)
    M[0][0]=2*(H[0]+H[1])
    M[0][1]=H[1]
    for j in range(1,n-2):
        M[j][j-1]=H[j]
        M[j][j]=2*(H[j]+H[j+1])
        M[j][j+1]=H[j+1]
    M[n-2][n-2]=2*(H[n-3]+H[n-2])
    M[n-2][n-3]=H[n-2]
    G=[]
    for q in range(n-1):
        g=-(6/H[q])*(F[q+1]-F[q])+(6/H[q+1])*(F[q+2]-F[q+1])
        G.append(g)
        
    return([M,G])
    
    
    
def zspln(M,G):
    '''input is the M an G as defined above,
    Matrix and Vector.
    The function outputs the Z values 
    where the first and last are both zero
    '''
    Z=[0]+list(np.linalg.inv(np.array(M)).dot(np.array(G)))+[0]
    return(Z)
    
    
 
    
def cspln(X,F,Z):
    '''input is the node list,
    f values at those nodes,
    and Z values at those nodes.
    outputs a matrix of coefficients for each
    spline used to interpolate the function.
    '''
    H=[]
    A=[]
    B=[]
    C=[]
    D=F.copy()
    for k in range(len(X)-1):
        H.append(X[k+1]-X[k])
    for i in range(len(X)-1):
        A.append((Z[i+1]-Z[i])/(6*H[i]))
        B.append(Z[i]/2)
        C.append((F[i+1]-F[i])/(H[i])-(H[i]*(Z[i+1]+2*Z[i]))/(6))
    Cmat=trans([A,B,C,D])
    return(Cmat)

def cube_array(COEF):
    coef_flip=[np.flip(v) for v in COEF]
    poly_list=[poly(u) for u in coef_flip]
    return(np.array(poly_list))

