#Code by GVV Sharma
#December 7, 2019
#Revised July 15, 2020
#released under GNU GPL
#Functions related to line
import numpy as np
import numpy.linalg as LA
import mpmath as mp
#from line.params import *
from params import *


def dir_vec(A,B):
  return B-A

def norm_vec(A,B):
  return omat@dir_vec(A,B)
  #return np.matmul(omat, dir_vec(A,B))

def ang_vec(m1,m2):
    return np.arccos(float((m1.T@m2)/(np.linalg.norm(m1)*np.linalg.norm(m2))))
    #return mp.acos(float((m1.T@m2)/(np.linalg.norm(m1)*np.linalg.norm(m2))))

#Generate line points
def line_gen(A,B):
  len =10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(0,1,len)
  for i in range(len):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB

#Generate line points
def line_gen_num(A,B,num):
  dim = A.shape[0]
  x_AB = np.zeros((dim,num))
  lam_1 = np.linspace(0,1,num)
  for i in range(num):
    temp1 = A + lam_1[i]*(B-A)
    x_AB[:,i]= temp1.T
  return x_AB
#Generating line in 2D using normal form
def line_norm(n,c,k1,k2):
    c = c/LA.norm(n)
    n = n/LA.norm(n)
    if c==0:
        A = np.zeros((2,1))
    elif np.array_equal(n, e1):
        A = np.array(([c, 0])).reshape(-1,1) 
    elif np.array_equal(n, e2):
        A = np.array(([0, c])).reshape(-1,1) 
    else:
        A = np.array(([c/n[0][0], 0])).reshape(-1,1) 
    m = omat@n
    return line_dir_pt(m,A,k1,k2)

#Normal to parametric
def param_norm(n,c):
    c = c/LA.norm(n)
    n = n/LA.norm(n)
    if c==0:
        A = np.zeros((2,1))
    elif np.array_equal(n, e1):
        A = np.array(([c, 0])).reshape(-1,1) 
    elif np.array_equal(n, e2):
        A = np.array(([0, c])).reshape(-1,1) 
    else:
        A = np.array(([c/n[0][0], 0])).reshape(-1,1) 
    m = omat@n
    return m,A

#Generating line using parametric form
def line_dir_pt(m,A,k1,k2):
  len = 10
  dim = A.shape[0]
  x_AB = np.zeros((dim,len))
  lam_1 = np.linspace(k1,k2,len)
  for i in range(len):
    temp1 = A + lam_1[i]*m
    x_AB[:,i]= temp1.T
  return x_AB

#Intersection of two lines
def line_isect(n1,c1,n2,c2):
  N=np.block([n1,n2]).T
  p = np.zeros((2,1))
  p[0] = c1
  p[1] = c2
  #Intersection
  P=np.linalg.solve(N,p)
  return P

#Intersection of two lines
def line_intersect(n1,A1,n2,A2):
  N=np.block([n1,n2]).T
  p = np.zeros((2,1))
  p[0] = n1.T@A1
  p[1] = n2.T@A2
  #Intersection
  P=np.linalg.solve(N,p)
  return P


#Foot of the perpendicular
def perp_foot(n,c,P):
  m = omat@n
  N=np.block([m,n])
  p = np.zeros((2,1))
  p[0][0] = m.T@P
  p[1][0] = c
  #Intersection
  x_0=np.linalg.solve(N.T,p)
  return x_0

#Rotation matrix
def rotmat(theta):
    #c = float(mp.cos(theta))
    #s = float(mp.sin(theta))
    c = np.cos(theta)
    s = np.sin(theta)
    P = np.array([[c,-s],[s,c]]) 
    return P
