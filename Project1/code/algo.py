import numpy as np
import matplotlib.pyplot as plt

f = lambda x: 100.0*np.exp(-10.0*x)
exact = lambda x: 1.0-(1-np.exp(-10))*x-np.exp(-10*x)
n = 10

#Matrix
b = np.zeros(n) #diagonal
a = np.zeros(n) #nedre diagonal
c = np.zeros(n) #øvre diagonal

# Multiplied with this vector
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
    b[i] = 2
    a[i] = -1
    c[i] = -1

exact_[0] = exact(x_0)
exact_[n] = exact(x_n)


for i in range(2,n):
    print(a[1] ==-1 and c[1] == -1)
    if a[1] !=-1 and c[1] == -1:
        g[i] = g[i] + g[i-1]/((i)/(i-1))
    else:
        b[i] = b[i] - a[i-1]*c[i-1]/b[i-1]
        g[i] = g[i] - a[i-1]*g[i-1]/b[i-1]


u[0] = 0 #boundry condition
u[n] = 0 #boundry condition
u[n-1] = g[n-1]/b[n-1] #boundry condition where u[n] = 0
i = n - 2
while i>0:
    if a[1] !=-1 and c[1] == -1:
        u[i] = (g[i] + u[i+1])/((i+1)/i) #special
    else:
        u[i] = (g[i] - c[i]*u[i+1])/b[i]
    i -= 1

"""
To test results do a quick plot or print out error
plt.plot(u, label = "computed")
plt.plot(exact_, linestyle = "dashed", label = "exact")
plt.legend()
plt.show()
"""

"""
for i in range(int(n/5)):
    print(f"x = {x0+i*h}:   ")
    #print(f"Calc: {u[i]}")
    #print(f"exact_: {exact_[i]}")
    print(f"Error: {abs(u[i]-exact_[i])}")
"""
