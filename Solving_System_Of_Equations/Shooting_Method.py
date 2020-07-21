# Shooting Method
# input function of th form y'' = p(x)y' + q(x)y + r(x)


import numpy as np
import matplotlib.pyplot as plt

def shootingMethod(p,q,r,a,b,alpha,beta,n):
    h = (b-a)/n     #sub interval size
    u = np.zeros([2,n+1])   #empty array
    v = np.zeros([2,n+1])   #empty array
    x = np.linspace(a,b,n+1) #cretes evenly spaced sequence

    u = solveSystem(fun1,fun2(p,q,r),a,alpha,0,b,n+1)  #runge-kutta-solve
    v = solveSystem(fun1,fun2(p,q,r),a,0,1,b,n+1)      #runge-kutta-solve
                        
    w = np.zeros([2,n+1])
    w[0,0] = alpha
    w[1,0] = (beta-u[0,n])/v[0,n]
    
    for i in range(1,n):
        w[0,i] = u[0,i] + w[1,0]*v[0,i]   
        w[1,i] = u[1,i] + w[1,0]*v[1,i]

    return(x,w)

#Runge-kutta-solving
# dy/dx = f1(x,y,z)  : differential equation 1
# dz/dx = f2(x,y,z)  : differential equation 2
# y(x(0)) = y0 : initial condition 1
# z(x(0)) = z0 : initial condition 2
# xf - variable final which define range
# n - number of points
def solveSystem(f1,f2,x0,y0,z0,xf,n):      
    h = (xf-x0)/(n-1)
    x = np.linspace(x0,xf,n)
    y = np.zeros([n,1])
    z = np.zeros([n,1])
    y[0] = y0
    z[0] = z0
    for i in range(1,n):
        k1 = h*f1(x[i-1],y[i-1],z[i-1])
        l1 = h*f2(x[i-1],y[i-1],z[i-1])
        k2 = h*f1(x[i-1]+h/2,y[i-1]+k1/2,z[i-1]+l1/2)
        l2 = h*f2(x[i-1]+h/2,y[i-1]+k1/2,z[i-1]+l1/2)
        k3 = h*f1(x[i-1]+h/2,y[i-1]+k2/2,z[i-1]+l2/2)
        l3 = h*f2(x[i-1]+h/2,y[i-1]+k2/2,z[i-1]+l2/2)
        k4 = h*f1(x[i-1]+h, y[i-1]+k3,z[i-1]+l3)
        l4 = h*f2(x[i-1]+h, y[i-1]+k3,z[i-1]+l3)
        
        y[i] = y[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6        
        z[i] = z[i-1] + (l1 + 2*l2 + 2*l3 + l4)/6
        
    return np.row_stack((np.transpose(y),np.transpose(z)))

def fun1(x, u, n):
    return n

def fun2(p, q, r):
    def fun22(x,v,n):
        return p(x)*n + q(x)*v + r(x)
    return fun22

def p(x):
    return 0

def q(x):
    return 3.34

def r(x):
    return 0

def main():
    a = 0   #stating point y(a) = alpha
    b = 4   #end point y(b) = beta
    alpha = 0.1   #boundry values
    beta = 0      #boundry values
    n = 40        #sub intervals
    x,w = shootingMethod(p,q,r,a,b,alpha,beta,n) # To use in general
    plt.plot(x,w[0],label='Y')
    #plt.plot(x,w[1],label='Y1')
    plt.xlabel("Length of tube in cm")
    plt.ylabel("Concentration in Molar(M)")
    plt.legend()
    plt.title("Solution with shooting method")
    plt.show()
    
main()