import numpy as np
import scipy.linalg

#########################################################
### y_dot = f(t, y);    y(t0) = y0
#########################################################
#h is the time step between two instants
#n is the dimension of y 

# y is an (n,m,N) shaped matrix, with y_i at time i in the ith column. This matrix includes y at all time instants. hence the M (all time instants)
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
        # k1 = k1/scipy.linalg.det(k1)
        k2=h*f(time[j]+h/2.0,y[:,:,j]+k1/2.0)
        # k2 = k2/scipy.linalg.det(k2)
        k3=h*f(time[j]+h/2.0,y[:,:,j]+k2/2.0)
        # k3 = k3/scipy.linalg.det(k3)
        k4=h*f(time[j]+h,y[:,:,j]+k3)
        # k4 = k4/scipy.linalg.det(k4)
        
        y[:,:,j+1]=y[:,:,j]+(k1+2*(k2+k3)+k4)/6.0
        # y[:,:,j+1]=y[:,:,j+1]/scipy.linalg.det(y[:,:,j+1])
        # print(k1/h)
        time[j+1]=time[j]+h
    x = y[:,:,-1]
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
    y=np.zeros((N,n))
    y[0,:] = y0


    for j in range(N-1):
        k1=h*f(time[j],y[j,:])
        k2=h*f(time[j]+h/2.0,y[j,:]+k1/2.0)
        k3=h*f(time[j]+h/2.0,y[j,:]+k2/2.0)
        k4=h*f(time[j]+h,y[j,:]+k3)
        # for i in range(n):
        #     k1[i]=h*f(time[j],y[:,j])[i]
        #     k2[i]=h*f(time[j]+h/2,y[:,j]+k1/2)[i]
        #     k3[i]=h*f(time[j]+h/2,y[:,j]+k2/2)[i]
        #     k4[i]=h*f(time[j]+h,y[:,j]+k3)[i]
        y[j+1,:]=y[j,:]+(k1+2*(k2+k3)+k4)/6.0
        time[j+1]=time[j]+h
    x = y[-1,:]
    # print(k1)
    return x, time



# Finds value of y for a given x using step size h 
# and initial value y0 at t0. 
# def RK4(f, t0, y0, x, h): 
#     # Count number of iterations using step size or 
#     # step height h 
#     n = (int)((x - t0)/h)  
#     # Iterate for number of iterations 
#     y = y0 
#     for i in range(1, n + 1): 
#         "Apply Runge Kutta Formulas to find next value of y"
#         k1 = h * f(t0, y) 
#         k2 = h * f(t0 + 0.5 * h, y + 0.5 * k1) 
#         k3 = h * f(t0 + 0.5 * h, y + 0.5 * k2) 
#         k4 = h * f(t0 + h, y + k3) 
  
#         # Update next value of y 
#         y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) 
  
#         # Update next value of x 
#         t0 = t0 + h 
#     return y 