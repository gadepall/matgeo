#Code by GVV Sharma
#December 7, 2019
#Revised November 12, 2023
#Revised July 15, 2020
#Revised August 1, 2024
#Revised August 14, 2024
#Revised August 17, 2024
#released under GNU GPL
#Functions related to conics

import numpy as np
import numpy.linalg as LA
from params import *

def conic_param(V,u,f):
        # Compute eigenvalues and eigenvectors
    lam,P = LA.eig(V)
    if lam[1]<=0 or lam[1] < lam[0]:
        lam = ref@lam
        P = ref@P
    e = np.sqrt(1-lam[0]/lam[1])
    p = P[:,0].reshape(-1,1)
    n = np.sqrt(np.abs(lam[1]))*p
    if e == 1:
        c = (LA.norm(u)**2-lam[1]*f)/(2*u.T@n)
        c = c[0][0]
        F = (c*e**2*n-u)/lam[1]
        #Standard parabola parameters
        eta = 2*u.T@p
        flen = -eta/lam[1]
        #Affine Parabola parameters
        cA = np.block([u+eta*p*0.5,V]).T
        cb = np.block([[-f],[0.5*eta*p-u]])
        O = LA.lstsq(cA,cb,rcond=None)[0]#vertex

    else:
        c = np.zeros((2,1)).flatten()#two directrices
        F = np.zeros((2,2))#two foci
        for i in range(2):
            disc = np.abs((e**2)*(u.T@n)**2-lam[1]*(e**2-1)*(LA.norm(u)**2-lam[1]*f))
            c[i] = (e*(u.T@n)+((-1)**i)*np.sqrt(disc))/(lam[1]*e*(e**2-1))
            F[:,i] = ((c[i]*(e**2)*n-u)/lam[1]).flatten()
        O = -LA.inv(V)@u
    return n,c,F,O,lam,P,e

#Standard parabola parameters
def parab_param(lam,P,u):
    p = P[:,0].reshape(-1,1)
    eta = 2*u.T@p
    flen = -eta/lam[1]
    return flen


#Standard ellipse/hyperbola parameters
def ellipse_param(V,u,f):
    lam,P = LA.eig(V)
    f0=u.T@LA.inv(V)@u-f
    ab =np.sqrt(np.abs(f0/lam))
    return ab.flatten()

#Generating points on a circle
def circ_gen(O,r):
	len = 50
	theta = np.linspace(0,2*np.pi,len)
	x_circ = np.zeros((2,len))
	x_circ[0,:] = r*np.cos(theta)
	x_circ[1,:] = r*np.sin(theta)
	x_circ = (x_circ + O)
	return x_circ

def circ_gen_num(O,r,num):
	theta = np.linspace(0,2*np.pi,num)
	x_circ = np.zeros((2,num))
	x_circ[0,:] = r*np.cos(theta)
	x_circ[1,:] = r*np.sin(theta)
	x_circ = (x_circ + O)
	return x_circ
#Generating points on an ellipse
def ellipse_gen(a,b):
	len = 50
	theta = np.linspace(0,2*np.pi,len)
	x_ellipse = np.zeros((2,len))
	x_ellipse[0,:] = a*np.cos(theta)
	x_ellipse[1,:] = b*np.sin(theta)
	return x_ellipse

def ellipse_gen_num(a,b,num):
	theta = np.linspace(0,2*np.pi,num)
	x_ellipse = np.zeros((2,num))
	x_ellipse[0,:] = a*np.cos(theta)
	x_ellipse[1,:] = b*np.sin(theta)
	return x_ellipse
#Generating points on a parabola
def parab_gen(y,a):
	x = y**2/a
	return x

#Generating points on a standard hyperbola 
def hyper_gen(y):
	x = np.sqrt(1+y**2)
	return x

#Points of contact for a pair of tangents
def contact(V,u,f,h):
    #intermediate
    gh = h.T@V@h+2*u.T@h+f 
    
    #matrix of tangents
    sigmat = (V@h+u)@(V@h+u).T-gh*V
    
    
    #Spectral decomposition
    D, P = LA.eig(sigmat)
    
    u1 = np.array(([np.sqrt(np.abs(D[1])),np.sqrt(np.abs(D[0]))]))
    u2 = np.array(([np.sqrt(np.abs(D[1])),-np.sqrt(np.abs(D[0]))]))
    
    u1 = u1.reshape(-1,1)
    u2 = u2.reshape(-1,1)
    
    #direction vectors
    m1 = P@u1
    m2 = P@u2
    #print(m1,m2)
    # Converting 1D array to a 2D numpy array of incompatible shape will cause error
    m1= np.reshape(m1, (2, 1))
    m2= np.reshape(m2, (2, 1))
    mu1 = -(m1.T@(V@h+u))/(m1.T@V@m1)
    mu2 = -(m2.T@(V@h+u))/(m2.T@V@m2)
    #print(mu1,mu2)
    x1 = h + mu1*m1
    x2 = h + mu2*m2
    return(x1,x2)
    
#Points of intersection for chords of a conic
def chord(V,u,f,m,h):
    #intermediate
    c = h.T@V@h+2*u.T@h+f 
    b = 2*m.T@(V@h+u)
    a = m.T@V@m
    k = np.roots((np.block([a,b,c])).flatten()).reshape(-1,1)
    Pmat = np.block([k,np.ones((2,1))]).T
    #print(Pmat,k)
    nmat = np.block([m,h])
    x = nmat@Pmat
    #x1 = h+k[0]*m
    #x2 = h+k[1]*m
    #return x1,x2
    return x
    
    
#Circle parameters
def circ_param(u,f):
    O = -u
    r = np.sqrt(LA.norm(u)**2-f)
    return O,r

#Circle tangent parameters
def circ_tang(n,u,r):
    Pmat = np.array(([1,-1],[-1,-1]))
    nmat = np.block([r*n,u])
    q = nmat@Pmat
    return q

#Conic tangent parameters
def conic_tangent(V,u,f,q):
    n = V@q+u
    c = -(u.T@q+f)
    #n = V@q+u
    #c = n.T@q
    return n,c

#Conic tangent contact
def conic_contact(V,u,f,n):
    _,_,_,_,lam,P,e = conic_param(V,u,f)
    p = P[:,0].reshape(-1,1)
    if e==1:
        #Parabola 
        r = (p.T@u).flatten()/(p.T@n).flatten()
        r = r[0]
        print(r)
        qA = np.block([u+r*n,V]).T
        qb = np.block([[-f],[r*n-u]])
        q = LA.lstsq(qA,qb,rcond=None)[0]#vertex
        print(qA,qb)
    else:
        #ellipse and hyperbola
        f0=u.T@LA.inv(V)@u-f
        r = np.sqrt(np.abs(f0/(n.T@LA.inv(V)@n)))
        Pmat = np.array(([1,-1],[-1,-1]))
        nmat = np.block([r*n,u])
        q = LA.inv(V)@nmat@Pmat
    return q
