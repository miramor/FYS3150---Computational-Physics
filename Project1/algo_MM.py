"""
Solving the one-dimensional Poisson equation with Dirchlet boundary conditions
Computing the solution with both a general algorithm and a specialised case, considering that the upper and lower diagonal are equal.
Algortihm includes forward and backward solution
"""

import numpy as np
import matplotlib.pyplot as plt
#
f = lambda x: 100.0*np.exp(-10.0*x)
#f = lambda x: 2*x+2
exact = lambda x: 1.0-(1-np.exp(-10))*x-np.exp(-10*x)
n = 10 #num of points

#Matrix
b = np.zeros(n) #diagonal
a = np.zeros(n) #lower diagonal
c = np.zeros(n) #upper diagonal

# Multiplied with this vector
u = np.zeros(n+1) # vector  x/v

#Results in this vector
g = np.zeros(n) # høyre side

#Stores exact solution to compare
exact_ = np.zeros(n+1)

x_0 = 0
x_n = 1
h = (x_n-x_0)/n # (x_n-x_0)/n
h_sq = h**2

#x_nplus1 = 0

for i in range(1,n):
    x_i = x_0+i*h
    exact_[i] = exact(x_i)
    g[i] = f(x_i)*h_sq
    a[i] = -1
    c[i] = -1
    if a[0] != c[0]: # if upper and lower diagonal not equal, b vector must be filled and cannot be precalculated
        b[i] = 2
    else:
        b[i] = (i+1.0)/i  #precalculate diagonal #2FLOPSx(n-1)

#exact values at boundary condition
exact_[0] = exact(x_0)
exact_[n] = exact(x_n)


for i in range(2,n):
    if  a[0] != c[0]:
         b[i] = b[i] - a[i-1]*c[i-1]/b[i-1] #3FLOPSx(n-2)
    g[i] = g[i] - a[i-1]*g[i-1]/b[i-1] #3FLOPsx(n-2)

#print(g)
#print(b)


u[0] = 0 #boundry condition
u[n] = 0 #boundry condition
u[n-1] = g[n-1]/b[n-1] #boundry condition where u[n] = 0 #1FLOP
i = n - 2
while i>0:
    u[i] = (g[i] - c[i]*u[i+1])/b[i] #3FLOPSx(n-2)
    i -= 1

# plt.plot(u, label = "computed")
# plt.plot(exact_, linestyle = "dashed", label = "exact")
# plt.legend()
# plt.show()

for i in range(1, n):
#print(exact_)
    error = np.log10(abs((u[i]-exact_[i])/exact_[i]))
# plt.plot((error))
# plt.show()
    print(error)

#print(f"Error: {abs(u[0:n:int(n/10)]-exact_[0:n:int(n/10)])}")

#Klart veldig mye bedre når c = a


"""for i in range(int(n/5)):
    print(f"x = {x0+i*h}:   ")
    #print(f"Calc: {u[i]}")
    #print(f"exact_: {exact_[i]}")
    print(f"Error: {abs(u[i]-exact_[i])}")"""
