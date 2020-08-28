import numpy as np
import matplotlib.pyplot as plt
#
f = lambda x: 100.0*np.exp(-10.0*x)
#f = lambda x: 2*x+2
exact = lambda x: 1.0-(1-np.exp(-10))*x-np.exp(-10*x)
n = 10 #num of points

#Matrix
b = np.zeros(n) #diagonal
a = np.zeros(n) #nedre diagonal
c = np.zeros(n) #øvre diagonal

# Multiplied with this vector
u = np.zeros(n+2) # vector  x/v

#Results in this vector
g = np.zeros(n) # høyre side

#Stores exact solution to compare
exact_ = np.zeros(n+2)

h = 1.0/n # (x_n-x_0)/n
h_sq = h**2

x0 = 0
#x_nplus1 = 0

for i in range(n): # 0, ______ , 0
    x_i = x0+(i+1)*h
    exact_[i+1] = exact(x_i)
    g[i] = f(x_i)*h_sq
    b[i] = 2
    a[i] = -1
    c[i] = -1

exact_[0] = exact(0)
exact_[n+1] = exact(x0+n*h)

print(x0+n*h)


for i in range(1,n):
    if a[0] != c[0]:
        #b[i] = -(i+1.0)/i
        b[i] = b[i] - 1/b[i-1]
    else:
        b[i] = b[i] - a[i-1]*c[i-1]/b[i-1]
    g[i] = g[i] - a[i-1]*g[i-1]/b[i-1]


u[n+1] = 0 #boundry condition where u[n+1] = 0
u[n] = g[n-1]/b[n-1] # second to last element in u is connected to last element in g, since g just contains the elements between the boundry condintions
i = n
while i>0:
    u[i] = (g[i-1] - c[i-1]*u[i+1])/b[i-1]
    i -= 1


#error = abs(u-exact_)/exact_
plt.plot(u, label = "Calculated")
plt.plot(exact_, linestyle = "dashed", label = "exact")
#plt.plot(abs(u[0:n+1:50]-exact_[0:n+1:50])/exact_[0:n+1:50])
plt.legend()
#plt.show()

#print(u[0])
#print(exact_[0])
#print(f"Error: {abs(u[0:n+1:50]-exact_[0:n+1:50])/exact_[0:n+1:50]}")

#Klart veldig mye bedre når c = a


"""for i in range(int(n/5)):
    print(f"x = {x0+i*h}:   ")
    #print(f"Calc: {u[i]}")
    #print(f"exact_: {exact_[i]}")
    print(f"Error: {abs(u[i]-exact_[i])}")"""
