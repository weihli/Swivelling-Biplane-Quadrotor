import numpy as np
import scipy.linalg

#########################################################
### y_dot = f(t, y);    y(t0) = y0
#########################################################
#h is the time step between two instants
#n is the dimension of y 

# y is an (n,N) shaped matrix, with y_i at time i in the ith column. This matrix includes y at all time instants. hence the M (all time instants)
def RK4_matrix(f, y0, t0, tf, h, n, m):        #h is the time step between two instants #n rows, m columns
    N=int((tf-t0)/h)
    k1=np.zeros((n,m))
    k2=np.zeros((n,m))
    k3=np.zeros((n,m))
    k4=np.zeros((n,m))
    time= np.linspace(t0, tf, N)

    y=np.zeros((n,m,N))
    y[:,:, 0] = y0


    for j in range(N-1):
        k1=h*f(time[j],y[:,:,j])
        k2=h*f(time[j]+h/2,y[:,:,j]+k1/2)
        k3=h*f(time[j]+h/2,y[:,:,j]+k2/2)
        k4=h*f(time[j]+h,y[:,:,j]+k3)
        
        y[:,:,j+1]=y[:,:,j]+(k1+2*(k2+k3)+k4)/6
        # y[:,:,j+1]=y[:,:,j+1]/scipy.linalg.det(y[:,:,j+1])
        time[j+1]=time[j]+h
    x = y[:,:,-1]
    # print(k1)
    return x, time


def RK4(f, y0, t0, tf, h, n):        #h is the time step between two instants
    N=int((tf-t0)/h)
    k1=np.array(np.zeros(n))
    k2=np.array(np.zeros(n))
    k3=np.array(np.zeros(n))
    k4=np.array(np.zeros(n))
    time= np.linspace(t0, tf, N)
    # time=np.array(np.zeros(N))
    # time[0]=t0
    y=np.array(np.zeros((n,N)))
    y[:, 0] = y0


    for j in range(N-1):
        k1=h*f(time[j],y[:,j])
        k2=h*f(time[j]+h/2,y[:,j]+k1/2)
        k3=h*f(time[j]+h/2,y[:,j]+k2/2)
        k4=h*f(time[j]+h,y[:,j]+k3)
        # for i in range(n):
        #     k1[i]=h*f(time[j],y[:,j])[i]
        #     k2[i]=h*f(time[j]+h/2,y[:,j]+k1/2)[i]
        #     k3[i]=h*f(time[j]+h/2,y[:,j]+k2/2)[i]
        #     k4[i]=h*f(time[j]+h,y[:,j]+k3)[i]
        y[:,j+1]=y[:,j]+(k1+2*(k2+k3)+k4)/6
        time[j+1]=time[j]+h
    x = y[:,-1]
    # print(k1)
    return x, time
