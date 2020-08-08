import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
#from RK4_matrix import RK4 as RK4

m = 0.8
g = 9.81
l = 0.42
T0 = m*g

end = 20
N = end*1000 + 1
time = np.linspace(0.0, end, N)
h = time[1] - time[0]
#states
R = np.zeros((N,3,3))
omg = np.zeros((N,3))
delt = np.zeros(N)
delt_d = np.zeros(N)

def hat(v):
	A = np.zeros((3,3))
	A[0,1] = -v[2]
	A[0,2] = v[1]
	A[1,2] = -v[0]
	A[1,0] = v[2]
	A[2,0] = -v[1]
	A[2,1] = v[0]
	return A

def vee(A):
	v = np.zeros(3)
	v[0] = A[2,1]
	v[1] = A[0,2]
	v[2] = A[1,0]
	return v

def axang2rotm(ax, ang):
    K = hat(ax/np.linalg.norm(ax))
    return np.eye(3) + np.sin(ang)*K + (1-np.cos(ang))*np.dot(K,K)

#initialize states
R[0,:,:] = axang2rotm(np.array([0.5, 0.0, 1.0]), np.pi)
#print(scipy.linalg.det(R[0,:,:]))
omg[0,:] = np.array([0.0,0.0,np.pi/2.0])
delt[0] = 0.0
delt_d[0] = 0.0

des_ax = np.array([1.0,1.0,1.0])
des_ax = des_ax/scipy.linalg.norm(des_ax)
freq = 2.0*np.pi
amp = 20.0*np.pi/180.0
axang_des = np.zeros((N,4))
omgdes = np.zeros((N,3))
omgdes_d = np.zeros((N,3))
omgdes_dd = np.zeros((N,3))
omgdes_ddd = np.zeros((N,3))


delt_des = np.zeros(N)
delt_desd = np.zeros(N)
delt_desdd = np.zeros(N)

Jxx = 0.01111
Jyy = 0.0136
Jzz = 0.02275

J = 2.0*np.array([[Jxx, 0.0, 0.0],
                  [0.0, Jyy, 0.0],
                  [0.0, 0.0, Jzz]])

J_inv = 0.5*np.array([[1/Jxx, 0.0, 0.0],
                      [0.0, 1/Jyy, 0.0],
                      [0.0, 0.0, 1/Jzz]])
    
k_R = 0.5*((2.0*np.pi)**2)*J
k_omg = 0.5*4.0*np.pi*J
eta_n = 8*np.pi
zeta = 1.0

D = 16*np.pi*np.eye(3)
K = (eta_n**2)*np.eye(3)

M = np.zeros((N,3))
M_d = np.zeros((N,3))
M_dd = np.zeros((N,3))
 
Mdes = np.zeros((N,3))
Mdes_d = np.zeros((N,3))
Mdes_dd = np.zeros((N,3))

Rdes = np.zeros((N,3,3))
Re = np.zeros((N,3,3))

omg_d = np.zeros((N,3))
omg_dd = np.zeros((N,3))
omg_ddd = np.zeros((N,3)) 

e_omg = np.zeros((N,3))
e_omgd = np.zeros((N,3))
e_omgdd = np.zeros((N,3))

phi = np.zeros(N)
th = np.zeros(N)
psi = np.zeros(N)

phi_d = np.zeros(N)
th_d = np.zeros(N)
psi_d = np.zeros(N)


A = np.array([[0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
              [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
              [-K[0,0], -K[0,1], -K[0,2], -D[0,0], -D[0,1], -D[0,2]],
              [-K[1,0], -K[1,1], -K[1,2], -D[1,0], -D[1,1], -D[1,2]],
              [-K[2,0], -K[2,1], -K[2,2], -D[2,0], -D[2,1], -D[2,2]]])

B = np.array([[0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [0.0, 0.0, 0.0],
              [1.0, 0.0, 0.0],
              [0.0, 1.0, 0.0],
              [0.0, 0.0, 1.0]])

M_arr = np.array([M[0,0], M[0,1], M[0,2], M_d[0,0], M_d[0,1], M_d[0,2]])
u1 = np.zeros(4) 
u1d = 0.0

def sw_biplane_model(R, omg, delt, delt_d, u, h):
    #R_del = np.array([[1.0, 0.0, 0.0],
    #                  [0.0, np.cos(delt), -np.sin(delt)],
    #                  [0.0, np.sin(delt), np.cos(delt)]])
    Ixx = 0.01111
    Iyy = 0.0136
    Izz = 0.02275
    I = 2.0*np.array([[Ixx, 0.0, 0.0],
                      [0.0, Iyy*np.cos(delt)**2 + Izz*np.sin(delt)**2, 0.0 ],
                      [0.0, 0.0, Izz*np.cos(delt)**2 + Iyy*np.sin(delt)**2]])
    I_inv = np.array([[1/I[0,0], 0.0, 0.0],
                      [0.0, 1/I[1,1], 0.0],
                      [0.0, 0.0, 1/I[2,2]]])
    I_dot = delt_d*(Iyy-Izz)*np.sin(delt)*np.array([[0.0,0.0,0.0],
                      [0.0,1.0,0.0],
                      [0.0,0.0,-1.0]])
    l = 0.42
    tau_m = u[0]
    tau_d = u[1]
    T_d = u[2]
    T_m = u[3]
    M = 2.0*np.array([tau_m, l*T_d*np.cos(delt),-l*T_m*np.sin(delt)])
    
    k1 = h*np.dot(I_inv, (-np.dot(I_dot, omg) - np.dot(hat(omg),np.dot(I,omg)) + M))
    k2 = h*np.dot(I_inv, (np.dot(-I_dot, (omg+k1/2.0)) - np.dot(hat(omg+k1/2.0),np.dot(I,(omg+k1/2.0))) + M))
    k3 = h*np.dot(I_inv, (np.dot(-I_dot, (omg+k2/2.0)) - np.dot(hat(omg+k2/2.0),np.dot(I,(omg+k2/2.0))) + M))
    k4 = h*np.dot(I_inv, (np.dot(-I_dot, (omg+k3)) - np.dot(hat(omg+k3),np.dot(I,(omg+k3))) + M))
    omg1 = omg + (k1+2.0*(k2+k3)+k4)/6.0
    delt_dd = (2.0*tau_d - np.sin(2*delt)*(Iyy-Izz)*(omg[1]**2 - omg[2]**2))/(2.0*Ixx)
    delt_d1 = delt_d + h*delt_dd
    #omg_del_1 = np.array([delt_d1, 0.0, 0.0])
    if omg1.all() == 0.0:
        R1 = R
    else:
        R1 = np.dot(R, axang2rotm(omg1/np.linalg.norm(omg1), np.linalg.norm(omg1)*h))
    #if omg_del_1.all() ==0:
    #    R_del_1 = R_del
    #else:
    #    R_del_1 = np.dot(R_del, axang2rotm(np.array([1.0,0.0,0.0]),np.linalg.norm(omg_del_1)*h))
    #del_1 = np.arctan(R_del_1[2,1]/R_del_1[1,1]) 
    del_1 = delt + delt_d1*h
    return R1, omg1, del_1, delt_d1

for i in range(N):
    t = time[i]
    ang = freq*t
    axang_des[i,:] = np.array([des_ax[0], des_ax[1], des_ax[2], amp*np.sin(ang)])
    omgdes[i,:] = amp*freq*np.cos(ang)*des_ax
    omgdes_d[i,:] = -amp*(freq**2)*np.sin(ang)*des_ax
    omgdes_dd[i,:] = -amp*(freq**3)*np.cos(ang)*des_ax
    omgdes_ddd[i,:] = amp*(freq**4)*np.sin(ang)*des_ax
    
    Rdes[i,:,:] = axang2rotm(des_ax, amp*np.sin(ang))
    Re[i,:,:] = np.dot(Rdes[i,:,:].transpose(), R[i,:,:])
    #print(Rdes[i,:,:])
    omg_d[i,:] = np.dot(J_inv, M[i,:]-np.dot(hat(omg[i,:]), np.dot(J, omg[i,:])))
    
    a1 = np.dot(hat(omg_d[i,:]), np.dot(J,omg[i,:]))
    a2 = np.dot(hat(omg[i,:]), np.dot(J,omg_d[i,:]))
    omg_dd[i,:] = np.dot(J_inv, M_d[i,:] - a1 - a2)
    
    e_omg[i,:] = omg[i,:] - np.dot(Re[i,:,:].transpose(), omgdes[i,:])
    e_omgd[i,:] = omg_d[i,:] - np.dot(Re[i,:,:].transpose(), omgdes_d[i,:]) + np.dot(hat(e_omg[i,:]), np.dot(Re[i,:].transpose(),omgdes[i,:]))
    
    A1 = np.dot(hat(e_omg[i,:]),np.dot(Re[i,:,:].transpose(),omgdes_d[i,:]))
    A2 = np.dot(Re[i,:,:].transpose(), omgdes_dd[i,:])
    A3 = np.dot(hat(e_omgd[i,:]),np.dot(Re[i,:,:].transpose(),omgdes[i,:]))
    A4 = np.dot(hat(e_omg[i,:]),np.dot(hat(e_omg[i,:]),np.dot(Re[i,:,:].transpose(),omgdes[i,:])))
    A5 = np.dot(hat(e_omg[i,:]),np.dot(Re[i,:,:].transpose(),omgdes_d[i,:]))
    e_omgdd[i,:] = omg_dd[i,:] + A1 - A2 + A3 - A4 + A5
    
    P = np.eye(3)
    eR = vee((np.dot(P, Re[i,:,:]) - np.dot(Re[i,:,:].transpose(), P)))/2.0
    eR_d = vee(np.dot(Re[i,:,:],hat(e_omg[i,:])) + np.dot(hat(e_omg[i,:]),Re[i,:,:].transpose()))/2.0
    
    B1 = np.dot(Re[i,:,:], np.dot(hat(e_omg[i,:]), hat(e_omg[i,:])))
    B2 = np.dot(hat(e_omgd[i,:]), Re[i,:,:].transpose())
    B3 = np.dot(np.dot(hat(e_omg[i,:]), hat(e_omg[i,:])), Re[i,:,:].transpose())
    eR_dd = vee(B1+B2-B3)/2.0
    
    C1 = np.dot(J, omgdes_d[i,:])
    C2 = np.dot(hat(omgdes[i,:]), np.dot(J, omgdes[i,:]))
    Mdes[i,:] = -np.dot(k_R,eR) - np.dot(k_omg,e_omg[i,:]) + np.dot(Re[i,:,:].transpose(), C1+C2)
    
    D1 = np.dot(hat(e_omg[i,:]),Re[i,:,:].transpose())
    D2 = np.dot(J, omgdes_d[i,:]) + np.dot(hat(omgdes[i,:]), np.dot(J, omgdes[i,:]))
    D3 = np.dot(J, omgdes_dd[i,:])
    D4 = np.dot(hat(omgdes_d[i,:]),np.dot(J,omgdes[i,:]))
    D5 = np.dot(hat(omgdes[i,:]),np.dot(J,omgdes_d[i,:]))
    Mdes_d[i,:] = -np.dot(k_R,eR_d) - np.dot(k_omg,e_omgd[i,:]) - np.dot(D1, D2) + np.dot(Re[i,:,:].transpose(),(D3 + D4 + D5))
    
    E1 = np.dot(hat(e_omg[i,:]),np.dot(hat(e_omg[i,:]),Re[i,:,:].transpose()))
    E2 = np.dot(hat(e_omgd[i,:]),Re[i,:,:].transpose())
    E3 = np.dot(J,omgdes_d[i,:]) + np.dot(hat(omgdes[i,:]), np.dot(J, omgdes[i,:]))
    E4 = np.dot(hat(e_omg[i,:]),Re[i,:,:].transpose())
    E5 = np.dot(J, omgdes_dd[i,:])
    E6 = np.dot(hat(omgdes_d[i,:]), np.dot(J, omgdes[i,:]))
    E7 = np.dot(hat(omgdes[i,:]),np.dot(J,omgdes_d[i,:]))
    E8 = np.dot(J,omgdes_ddd[i,:])
    E9 = np.dot(hat(omgdes_dd[i,:]), np.dot(J, omgdes[i,:]))
    E10 = 2.0*np.dot(hat(omgdes_d[i,:]), np.dot(J, omgdes_d[i,:]))
    E11 = np.dot(hat(omgdes[i,:]),np.dot(J, omgdes_dd[i,:]))
    Mdes_dd[i,:] = -np.dot(k_R,eR_dd) - np.dot(k_omg,e_omgdd[i,:]) + np.dot((E1-E2),E3)- np.dot(E4, (E5 + E6 + E7)) + np.dot(Re[i,:,:].transpose(), (E8 + E9 + E10 + E11))
    
    u = Mdes_dd[i,:] + np.dot(D, Mdes_d[i,:]) + np.dot(K, Mdes[i,:])

    M_arr = M_arr + (np.dot(A, M_arr) + np.dot(B, u))*h
    if i<N-1:
        M[i+1,:] = M_arr[:3]
        M_d[i+1,:] = M_arr[3:]
    
    M[i,2] = -2.0*l*T0*np.tan(delt[i])
    M_d[i,2] = -2.0*l*T0*((1.0/np.cos(delt[i]))**2)*delt_d[i]
    #M[i,:] = M_arr[:3]
    #M_d[i,:] = M_arr[3:]
    
    delt_des[i] = -np.arctan(Mdes[i,2]/(2.0*l*T0))
    delt_desd[i] = -2.0*l*T0*Mdes_d[i,2]/(Mdes[i,2]**2 + (2*l*T0)**2)
    delt_desdd[i] = -2.0*l*T0*Mdes_dd[i,2]/(Mdes[i,2]**2 + (2*l*T0)**2) + 4.0*l*T0*Mdes[i,2]*(Mdes_d[i,2]**2)/((Mdes[i,2]**2 + (2.0*l*T0)**2)**2)
    
    v1 = delt_desdd[i] - 2.0*2.0*np.pi*2.0*(delt_d[i] - delt_desd[i]) - ((2*np.pi*2)**2)*(delt[i] - delt_des[i])
    #print(v1*h*h*180.0/np.pi)
    c1 = np.sin(2.0*delt[i])*(Jyy-Jzz)*((omg[i,1])**2 - (omg[i,2])**2)/(2.0*Jxx)
    u1[1] = Jxx*(c1+v1)   ### tau_delta
    u1[0] = M_arr[0]/2.0
    u1[2] = M_arr[1]/(2.0*l*np.cos(delt[i]))
    u1[3] = T0/np.cos(delt[i])
    
    u1d = u1d + h*(u1-u1d)/(0.015)
    #u1 = u - np.dot(D, M_d[i,:]) - np.dot(K, M[i,:])
    #nu_z =  -u1[2]*(np.cos(delt[i])**2)/(2.0*l*T0) - 2.0*(delt_d[i]**2)*np.tan(delt[i])
    #tau_delta = Jxx*v1 + 0.5*(Jyy-Jzz)*np.sin(2*delt[i])*((omg[i,1]**2)-(omg[i,2]**2))
    #print(i, "      ", omg[i,:]*180.0/np.pi)
    if i<N-1:
        R[i+1,:,:] = sw_biplane_model(R[i,:,:], omg[i,:], delt[i], delt_d[i], u1d, h)[0]
        omg[i+1,:] = sw_biplane_model(R[i,:,:], omg[i,:], delt[i], delt_d[i], u1d, h)[1]
        delt[i+1] = sw_biplane_model(R[i,:,:], omg[i,:], delt[i], delt_d[i], u1d, h)[2]
        delt_d[i+1] = sw_biplane_model(R[i,:,:], omg[i,:], delt[i], delt_d[i], u1d, h)[3]
    
    #print(delt[i])
    phi[i] = np.arcsin(R[i,2,1])
    th[i]  = np.arctan(-R[i,2,0]/R[i,2,2])
    psi[i] = np.arctan(-R[i,0,1]/R[i,1,1])
    #print(R[i,1,1])
    
    phi_d[i] = np.arcsin(Rdes[i,2,1])
    th_d[i]  = np.arctan(-Rdes[i,2,0]/Rdes[i,2,2])
    psi_d[i] = np.arctan(-Rdes[i,0,1]/Rdes[i,1,1])
    
fig, ax = plt.subplots(2,2)
ax[0,0].plot(time, 180.0*(1/np.pi)*phi,'r')
ax[0,0].plot(time, 180.0*(1/np.pi)*phi_d,'b')
ax[0,0].set_xlabel('Time')
ax[0,0].set_ylabel('Phi')
ax[0,0].legend(['phi', 'des phi'], fontsize='xx-small')

ax[0,1].plot(time, 180.0*(1/np.pi)*th,'r')
ax[0,1].plot(time, 180.0*(1/np.pi)*th_d,'b')
ax[0,1].set_xlabel('Time')
ax[0,1].set_ylabel('Theta')
ax[0,1].legend(['theta', 'des theta'], fontsize='xx-small')

ax[1,0].plot(time, 180.0*(1/np.pi)*psi,'r')
ax[1,0].plot(time, 180.0*(1/np.pi)*psi_d,'b')
ax[1,0].set_xlabel('Time')
ax[1,0].set_ylabel('Psi')
ax[1,0].legend(['psi', 'des psi'], fontsize='xx-small')

ax[1,1].plot(time, 180.0*(1/np.pi)*delt,'r')
ax[1,1].plot(time, 180.0*(1/np.pi)*delt_des,'b')
ax[1,1].set_xlabel('Time')
ax[1,1].set_ylabel('delta')
ax[1,1].legend(['delta', 'delta_des'], fontsize='xx-small')

plt.show()