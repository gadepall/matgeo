#Code by GVV Sharma
#December 7, 2019
#Revised July 15, 2020
#Revised October 21, 2023
#released under GNU GPL
#Functions related to plots

import numpy as np
import numpy.linalg as LA
from line.funcs import *
from params import *
import matplotlib.pyplot as plt


#Labeling points
def label_pts(G_v,vert_labels): 
    for i, txt in enumerate(vert_labels):
        plt.annotate(txt, # this is the text
                     (G_v[0,i], G_v[1,i]), # this is the point to label
                     textcoords="offset points", # how to position the text
                     xytext=(0,10), # distance from text to points (x,y)
                     ha='center') # horizontal alignment can be left, right or center


