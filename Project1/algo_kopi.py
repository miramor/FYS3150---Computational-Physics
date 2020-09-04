import numpy as np
import matplotlib.pyplot as plt
#
f = lambda x: 100.0*np.exp(-10.0*x)
#f = lambda x: 2*x+2
exact = lambda x: 1.0-(1-np.exp(-10))*x-np.exp(-10*x)
n = 100 #num of points

#Matrix
b = np.zeros(n) #diagonal
a = np.zeros(n) #nedre diagonal
c = np.zeros(n) #øvre diagonal

# Multiplied with this vector
u = np.zeros(n) # vector  x/v

#Results in this vector
g = np.zeros(n) # høyre side

#Stores exact solution to compare
exact_ = np.zeros(n)

h = 1.0/n # (x_n-x_0)/n
h_sq = h**2

x0 = 0
#x_nplus1 = 0

for i in range(n):
    x_i = x0+i*h
    exact_[i] = exact(x_i)
    g[i] = f(x_i)*h_sq
    b[i] = 2
    a[i] = -1
    c[i] = -1

for i in range(1,n):
    if a[0] == c[0]:
        #b[i] = -(i+1.0)/i
        b[i] = b[i] - 1/b[i-1]
    else:
        b[i] = b[i] - a[i-1]*c[i-1]/b[i-1]
    g[i] = g[i] - a[i-1]*g[i-1]/b[i-1]

#print(g)
#print(b)


u[n-1] = g[n-1]/b[n-1] #boundry condition where u[n] = 0
i = n - 2
while i>=0:
    u[i] = (g[i] - c[i]*u[i+1])/b[i]
    i -= 1

plt.plot(u, label = "computed")
plt.plot(exact_, label = "exact")
plt.legend()
plt.show()


#print(f"Error: {abs(u[0:n:50]-exact_[0:n:50])}")

#Klart veldig mye bedre når c = a


"""for i in range(int(n/5)):
    print(f"x = {x0+i*h}:   ")
    #print(f"Calc: {u[i]}")
    #print(f"exact_: {exact_[i]}")
    print(f"Error: {abs(u[i]-exact_[i])}")"""
