import numpy as np
import  matplotlib.pyplot as plt
import scipy.linalg
from RK4_matrix import RK4 as RK4
from RK4_matrix import RK4_matrix as RK4_matrix
# from RK4 import RK4

N = 201
end = 20.0
time = np.linspace(0.0, end, N)
R_d = np.zeros((N,3,3))
R = np.zeros((N,3,3))

delta = np.zeros(N)
delta[0] = (np.pi/180.0)*10.0

Jxx = 0.01111
Jyy = 0.0136
Jzz = 0.02275

phi = np.zeros(N)
th = np.zeros(N)
psi = np.zeros(N)

phi_d = np.zeros(N)
th_d = np.zeros(N)
psi_d = np.zeros(N)
R_d_dot = np.zeros((N,3,3))
R_dot = np.zeros((N,3,3))
R_e = np.zeros((N,3,3))
R_e[0,:,:] = np.eye(3)
R_e_dot = np.zeros((N,3,3))
e_R = np.zeros((N,3,3))
P = np.array([[1.0,0.0,0.0],
			  [0.0,1.0,0.0],
			  [0.0,0.0,1.0]])

k_R = 0.1
k_w = 0.1

Moment = np.zeros((N,3))
Moment_dot = np.zeros((N,3))

M_d = np.zeros((N,3))

M_e = np.zeros((N,3))
M_e[0,:] = np.array([0.1, 0.2, 0.3])

M_e_dot = np.zeros((N,3))
M_e_dot[0,:] = np.array([0.01,0.2,0.2])

X = np.zeros((N,6))
X[0,:] = np.concatenate((M_e[0,:], M_e_dot[0,:]))

xi    = 2.0*np.array([0.01,0.01,0.01])
omega = np.array([0.01,0.01,0.01])

D = 2.0*np.array([[xi[0]*omega[0],0.0,0.0],
	              [0.0,xi[1]*omega[1],0.0],
	              [0.0,0.0,xi[2]*omega[2]]])

K = np.array([[omega[0]**2,0.0,0.0],
			  [0.0,omega[1]**2,0.0],
			  [0.0,0.0,omega[2]**2]])
M = np.zeros((6,6))
M[:3,3:] = np.eye(3)
M[3:,:3] = K
M[3:,3:] = D

for i in range(N-1):
	X[i+1,:] = np.dot(scipy.linalg.expm(time[i]*M), X[0])

for i in range(N):
	M_e[i,:] = X[i, :3]
	M_e_dot[i,:] = X[i,3:]

phi_d_dot = np.zeros(N)
th_d_dot = np.zeros(N)
psi_d_dot = np.zeros(N)


w = np.zeros((N, 3))
w[0,:] = (np.pi/180.0)*np.array([10.0,10.0,10.0])
w_d = np.zeros((N,3))
w_d_dot = np.zeros((N,3))
e_w = np.zeros((N,3))


def rotm(phi, th, psi):
	R = np.zeros((3,3))
	R[0,0] = np.cos(th)*np.cos(psi) - np.sin(phi)*np.sin(th)*np.sin(psi)
	R[0,1] = -np.cos(phi)*np.sin(psi)
	R[0,2] = np.sin(th)*np.cos(psi) + np.sin(phi)*np.cos(th)*np.sin(psi)
	R[1,0] = np.cos(th)*np.sin(psi) + np.sin(phi)*np.sin(th)*np.cos(psi)
	R[1,1] = np.cos(phi)*np.cos(psi)
	R[1,2] = np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(th)*np.cos(psi)
	R[2,0] = -np.cos(phi)*np.sin(th)
	R[2,1] = np.sin(phi)
	R[2,2] = np.cos(phi)*np.cos(th)
	return R

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


R[0,:,:] = rotm(0.0, 50.0*np.pi/180.0, np.pi)
R_d[0,:,:] = np.eye(3)
R_e[0,:,:] = np.dot(R_d[0,:,:].transpose(),R[0,:,:])

# print(R[0,:,:])

for i in range(N):

	J = 2.0*np.array([[Jxx, 0.0, 0.0],
				  [0.0, Jyy, 0.0 ],
				  [0.0, 0.0, Jzz*np.cos(delta[i])**2 + Jyy*np.sin(delta[i])**2]])


	# print(i)

	def rotm_int_w(t, R):
		return np.dot(R, hat(w[i,:]))

	def rotm_int_ew(t, R):
		return np.dot(R, hat(e_w[i,:]))

	phi_d[i] = (np.pi/180.0)*20.0*np.sin(2*np.pi*time[i]/1.0)
	th_d[i] = (np.pi/180.0)*20.0*np.sin(2*np.pi*time[i]/1.0)
	psi_d[i] = (np.pi/180.0)*20.0*np.sin(2*np.pi*time[i]/1.0)
	R_d[i,:,:] = rotm(phi_d[i], th_d[i], psi_d[i])
	phi_d_dot[i] = (np.pi/180.0)*(2*np.pi)*20.0*np.cos(2*np.pi*time[i]/1.0)
	th_d_dot[i] = (np.pi/180.0)*(2*np.pi)*20.0*np.cos(2*np.pi*time[i]/1.0)
	psi_d_dot[i] = (np.pi/180.0)*(2*np.pi)*20.0*np.cos(2*np.pi*time[i]/1.0)
	w_d[i,0] = phi_d_dot[i]*np.sin(th_d[i])*np.sin(psi_d[i]) + th_d_dot[i]*np.cos(psi_d[i])
	w_d[i,1] = phi_d_dot[i]*np.sin(th_d[i])*np.cos(psi_d[i]) - th_d_dot[i]*np.sin(psi_d[i])
	w_d[i,2] = phi_d_dot[i]*np.cos(th_d[i]) + psi_d_dot[i]

	e_w[0,:] = w[0,:]- np.dot(R_e[0,:,:].transpose(),w_d[0,:])

	# R_e[i,:,:] = np.dot(R_d[i,:,:].transpose(), R[i,:,:])
	# e_w[i,:] = w[i,:] - np.dot(R_e[i,:,:].transpose(), w_d[i,:])
	e_R[i,:,:] = 0.5*(np.dot(P,R_e[i,:,:]) - np.dot(R_e[i,:,:].transpose(),P))

	# print(M_e[i,:])


	def error_dyn(t, e_w):
		return np.dot(scipy.linalg.inv(J), -k_R*vee(e_R[i,:,:]) -k_w*e_w + M_e[i,:])


	if i<N-1:
		e_w[i+1,:] = RK4(error_dyn, e_w[i,:], time[i], time[i+1], 0.001, 3)[0]

		R_e[i+1,:,:] = RK4_matrix(rotm_int_ew, R_e[0,:,:], time[i], time[i+1], 0.001, 3, 3)[0]
		R_e[i+1,:,:] = R_e[i+1,:,:]/scipy.linalg.det(R_e[i+1,:,:])

		w_d_dot[i,:] = (w_d[i+1,:]-w_d[i,:])/(time[i+1]-time[i])

	R[i,:,:] = np.dot(R_d[i,:,:], R_e[i,:,:])
	w[i,:] = e_w[i,:] + np.dot(R_e[i,:,:].transpose(), w_d[i,:])

	var = np.dot(hat(e_w[i,:]), np.dot(R_e[i,:,:].transpose(), w_d[i,:])) - np.dot(R_e[i,:,:].transpose(), w_d_dot[i,:])
	M_d[i,:] = -k_R*vee(e_R[i,:,:]) - k_w*e_w[i,:] + np.cross(w[i,:], np.dot(J,w[i,:])) - np.dot(J, var)

	Moment[i,:] =  M_e[i,:] + M_d[i,:]

	if i<N-1:
		Moment_dot[i,:] = (Moment[i+1]-Moment[i])/(time[i+1]-time[i])

	# print(Moment_dot[i,:])

	phi[i] = np.arcsin(R[i,2,1])
	th[i]  = np.arccos(R[i,2,2]/np.cos(phi[i]))
	psi[i] = np.arccos(R[i,1,1]/np.cos(phi[i]))

	def delta_int(t, delta):
		return (Moment_dot[i,2]/Moment[i,2])*np.sin(2*delta)

	if i<N-1:
		delta[i+1] = RK4(delta_int, delta[i], time[i], time[i+1], 0.001, 1)[0]
	# print(delta[i])


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

ax[1,1].plot(time, 180.0*(1/np.pi)*delta,'r')
# ax[0,0].plot(time, 180.0*(1/np.pi)*psi_d,'r', linestyle='dashed')
ax[1,1].set_xlabel('Time')
ax[1,1].set_ylabel('delta')
ax[1,1].legend(['delta'], fontsize='xx-small')

plt.show()