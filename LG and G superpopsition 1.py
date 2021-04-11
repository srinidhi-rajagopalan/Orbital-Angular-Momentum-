import numpy #Import numerical python library 
import matplotlib.pyplot as plt #Import library for plotting 
import math
from scipy import special
w=100 #peak to peak diameter of the doughnut
x= numpy.linspace(-100,100,30) #Makes an array of x coordinate 
y= numpy.linspace(-100,100,30) #Makes an array of y coordinate 
shape = 30 #Defines the number of elements in your array
lamda=8.4 
l=1
n=1.84
d=10

#-------------------------LAGUERE GAUSSIAN-----------------------------------------
def LG(x,y): #The function describes the LG mode 
    return (2/numpy.pi**0.5) *(1/w**2) * ((x**2+y**2)**0.5)*numpy.exp(-(x**2+y**2)/w**2)#*numpy.exp (1/(1+(y/x)**2))
    #return (2**0.5)/w  *(x**2+y**2)**0.5 *numpy.exp(-(x**2+y**2)/w**2)
LG_array=[LG(elemx,elemy)for elemx in x for elemy in y] #Solves the function for x and y elements and it is of type'list'
LG_array=numpy.array(LG_array) #This converts the 'list' to 'array'. This is now a column of 900 elements
LG_array=LG_array.reshape (shape,shape) #This converts the column into a 30X30 matrix
plt.pcolormesh(x,y,LG_array,cmap='PuBu_r') #Function for plotting into color distribution
plt.colorbar() #This includes a colorbar to the right of your plot
plt.show()
#--------------------------GAUSSIAN---------------------------------------

def G(x,y):
    return numpy.exp(-(x**2+y**2)/w**2)
G_array=[G(elemx,elemy)for elemx in x for elemy in y] #Solves the function for x and y elements and it is of type'list'
G_array=numpy.array(G_array) #This converts the 'list' to 'array'. This is now a column of 900 elements
G_array=G_array.reshape (shape,shape) #This converts the column into a 30X30 matrix
plt.pcolormesh(x,y,G_array,cmap='PuBu_r') #Function for plotting into color distribution
plt.colorbar()
plt.show()
#---------------------------INTERFERENCE PDF --------------------------------------
'''
def teta_0(L):
    return 2* numpy.pi * L/(lamda * l) 
def teta(L):
    return (2*numpy.pi *n/l) - teta_0(L)
def dphi(L):
    return l*teta(L) - ((2*numpy.pi * L)/lamda)
def final(x,y,L):
    return 2*LG(x,y)*G(x,y) *(1+numpy.cos(dphi(L)))
L=11
print (dphi(L))
'''
#----------------------INTERFERENCE 2--------------------------------------------------
def dphi(x,y,L):
    return (2*numpy.pi *L)/lamda + (2*numpy.pi/lamda) * d * numpy.arctan2(y,x)/8
def final(x,y,L):
    return    2*LG(x,y)*G(x,y) *(1+numpy.cos(dphi(x,y,L)))
L=26
print (dphi(x,y,L))

final=[final(elemx,elemy,L)for elemx in x for elemy in y] #Solves the function for x and y elements and it is of type'list'
final=numpy.array(final) #This converts the 'list' to 'array'. This is now a column of 900 elements
final=final.reshape (shape,shape) #This converts the column into a 30X30 matrix
plt.pcolormesh(x,y,final,cmap='PuBu_r') #Function for plotting into color distribution
plt.show() #Function that seperates and displays the three different plots








   

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

