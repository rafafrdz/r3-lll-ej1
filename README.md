![](https://raw.githubusercontent.com/rafafrdz/r3-lll-ej1/master/img/R_1.png)

![](https://raw.githubusercontent.com/rafafrdz/r3-lll-ej1/master/img/R_2.png)

```python
# -*- coding: utf-8 -*-
import time
from numpy import *
from scipy.sparse.linalg import cgs,splu
from scipy.sparse import lil_matrix,identity
from scipy.linalg import lu_factor,lu_solve,cho_factor,cho_solve
from matplotlib.pyplot import *

def f0(x):
    return 0*x
def u0(x):
    y = exp(-(x-5)**2)
    return y

def ejercicio(x0,xf,N,T,Nt,mu,ua,ub,fuente,cond0):
    t1=time.time()
    N = int(N)
    Nt = int(Nt)
    x0 = float(x0)
    xf = float(xf)
    x = linspace(x0,xf,N+1)
    dx = (xf-x0)/float(N)
    dx2 = dx**2
    dt = T/float(Nt)
    #Configuramos M
    M = lil_matrix((N+1,N+1),dtype='float64')
    M.setdiag(2.0*ones(N+1),0)
    M.setdiag(-1.0*ones(N+1),1)
    M.setdiag(-1.0*ones(N+1),-1)
    M[0,0] = 0.0
    M[N,N] = 0.0
    M[0,1]= 0.0
    M[1,0]= 0.0
    M[N,N-1]=0.0
    M[N-1,N]= 0.0
    M=M.tocsc()
    #Configuramos la matriz Identidad
    Id = identity(N+1, dtype='float64', format='csc')
    #Damos forma a la matriz A
    A = Id + mu*dt/dx2*M
    #Configuramos el vector B
    b = fuente(x)
    #Imponemos la condicion inicial
    usol = cond0(x)
    #hace la descomposición LU completa de una matriz Sparse
    LU=splu(A)
    #plot(x,usol,'r')
    for i in range(Nt):
        b = dt*fuente(x) + usol
        b[0]=ua
        b[1]+= mu*dt*ua/dx2
        b[N-1] += mu*dt*ub/dx2
        b[N]=ub
        usol=LU.solve(b)
        #plot(x,usol,'g')
    tf=time.time()
    print "Tiempo de ejecucion:",tf-t1
    plot(x,usol,'b')
    show()

#ejercicio(0,20,200,20,400, 0.1, 0, 0, f0, u0)
#ejercicio(0,20,200,20,380, 0.1, 0, 0, f0, u0)
#ejercicio(0,20,200,20,100, 0.1, 0, 0, f0, u0)

```

