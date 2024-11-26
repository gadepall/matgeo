#Code by GVV Sharma
#December 7, 2019
#Revised July 15, 2020
#Revised May 15, 2024
#released under GNU GPL
#Functions related to triangle

import numpy as np
import numpy.linalg as LA
from line.funcs import *
#from line.params import *
from params import *

#Triangle sides
def tri_sides(A,B,C):
    a = LA.norm(B-C)
    b = LA.norm(C-A)
    c = LA.norm(A-B)
    return c,a,b

#Triangle Construction 
#a,b,c
def tri_vert(a,b,c):
  p = (a**2 + c**2-b**2 )/(2*a)
  q = np.sqrt(c**2-p**2)
  A = np.array([p,q]) 
  B = np.array([0,0]) 
  C = np.array([a,0]) 
  return  A,B,C

#a,b+c,angB given
def tri_const(a,K,angB):
  c= 0.5*(K**2-a**2)/(K-a*np.cos(angB))
  C=np.array(([a,0])).reshape(-1,1)
  B=np.array(([0,0])).reshape(-1,1)
  A=np.array(([c*np.cos(angB),c*np.sin(angB)])).reshape(-1,1)
  return  A,B,C

#a+b+c,angB,angC given
def tri_const_sum(K,angB,angC):
  x=LA.solve(np.array(([1,1,1],[-1,np.cos(angC),np.cos(angB)],[0,np.sin(angC),-np.sin(angB)])),K*np.array(([1,0,0])).reshape(-1,1))
  print(x)
  C=np.array(([x[0][0],0])).reshape(-1,1)
  B=np.array(([0,0])).reshape(-1,1)
  A=np.array(([x[2][0]*np.cos(angB),x[2][0]*np.sin(angB)])).reshape(-1,1)
  return  A,B,C

#Triangle  mid points
def tri_mid_pt(A,B,C):
  D = (B+C)/2
  E = (C+A)/2
  F = (A+B)/2
  return  D,E,F


#Radius and centre of the circumcircle
#of triangle ABC
def ccircle(A,B,C):
    D,E,F = tri_mid_pt(A,B,C)
    m1 = dir_vec(A,B)
    m2 = dir_vec(A,C)
    O = line_intersect(m1,F,m2,E)
    r = LA.norm(A -O)
    return O,r

#Radius and centre of the incircle
#of triangle ABC
def icircle(A,B,C):
  k1 = 1
  k2 = 1
  p = np.zeros((2,1))
  t = norm_vec(B,C)
  n1 = t/LA.norm(t)
  t = norm_vec(C,A)
  n2 = t/LA.norm(t)
  t = norm_vec(A,B)
  n3 = t/LA.norm(t)
  p[0] = n1.T@B- k1*n2.T@C
  p[1] = n2.T@C- k2*n3.T@A
  #Intersection
  N1= n1 - k1 * n2
  N2= n2 - k1 * n3
  N=np.block([[N1.T],[N2.T]])
  I=LA.solve(N,p)
  r = n1.T@(I-B)
  return I,r



#Incircle points of contact
'''
def icontact(A,B,C):
    LA.inv(I_circ_mat)@tri_sides(A,B,C).(-1,1)
'''
