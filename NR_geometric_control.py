import numpy as np
import  matplotlib.pyplot as plt
import scipy.linalg
from RK4_matrix import RK4 as RK4
from RK4_matrix import RK4_matrix as RK4_matrix
# from RK4 import RK4

m = 0.8 #kg
g = 9.81
l = 0.42 #m

N = 201
end = 20.0
time = np.linspace(0.0, end, N)
R_d = np.zeros((N,3,3))
R = np.zeros((N,3,3))

delta = np.zeros(N)
delta[0] = (np.pi/180.0)*10.0

delta_dot = np.zeros(N)

delta_des = np.zeros(N)
delta_des[0] = (np.pi/180.0)*1.0

delta_des_dot = np.zeros(N)

Jxx = 0.01111
Jyy = 0.0136
Jzz = 0.0136#0.02275



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

k_R = 0.8
k_w = 0.8
Moment = np.zeros((N,3))
Moment_dot = np.zeros((N,3))
Moment_dot_dot = np.zeros((N,3))

M_d = np.zeros((N,3))

############ Moment error calculation begins###############

M_e = np.zeros((N,3))
M_e[0,:] = np.array([0.1, 0.2, 0.3])

M_e_dot = np.zeros((N,3))
M_e_dot[0,:] = np.array([0.01,0.2,0.2])

X = np.zeros((N,6))
X[0,:] = np.concatenate((M_e[0,:], M_e_dot[0,:]))

xi    = np.array([1.0,1.0,1.0])
omega =2*np.pi*np.array([1.0,1.0,1.0])

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
	X[i+1,:] = np.dot(scipy.linalg.expm(time[i+1]*M), X[0])

for i in range(N):
	M_e[i,:] = X[i, :3]
	M_e_dot[i,:] = X[i,3:]

########## Moment error calculation ends ###############

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

def rodrigues(omega, t):
	return np.eye(3) + hat(omega/scipy.linalg.norm(omega))*np.sin(t*scipy.linalg.norm(omega)) + np.dot(hat(omega/scipy.linalg.norm(omega)),hat(omega/scipy.linalg.norm(omega)))*(1-np.cos(t*scipy.linalg.norm(omega)))


R[0,:,:] = rotm(0.0, 50.0*np.pi/180.0, np.pi)
R_d[0,:,:] = np.eye(3)
R_e[0,:,:] = np.dot(R_d[0,:,:].transpose(),R[0,:,:])

M_d_dot = np.zeros((N, 3))
M_d_dot_dot = np.zeros((N,3))

T_m = np.zeros((N,3))

x = np.zeros((N,2))
x[0,:] = np.array([delta[0], delta_dot[0]])
y = np.zeros((N,2))
y[0,:] = np.array([delta_des[0], delta_des_dot[0]])
# print(scipy.linalg.norm(np.array([0.0,0.0,m*g])))
# print(R[0,:,:])

for i in range(N):

	J = 2.0*np.array([[Jxx, 0.0, 0.0],
				      [0.0, Jyy, 0.0 ],
				      [0.0, 0.0, Jzz]])

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

	R_e[0,:,:] = np.dot(R_d[0,:,:].transpose(), R[0,:,:])
	# e_w[i,:] = w[i,:] - np.dot(R_e[i,:,:].transpose(), w_d[i,:])
	e_R[i,:,:] = 0.5*(np.dot(P,R_e[i,:,:]) - np.dot(R_e[i,:,:].transpose(),P))

	# print(e_R[i,:,:])

	# print(M_e[i,:])


	def error_dyn(t, e_w):
		return np.dot(scipy.linalg.inv(J), -k_R*vee(e_R[i,:,:]) -k_w*e_w + M_e[i,:])


	if i<N-1:
		e_w[i+1,:] = RK4(error_dyn, e_w[i,:], time[i], time[i+1], 0.001, 3)[0]

		# R_e[i+1,:,:] = RK4_matrix(rotm_int_ew, R_e[0,:,:], time[i], time[i+1], 0.001, 3, 3)[0]
		R_e[i+1,:,:] = np.dot(R_e[i,:,:], rodrigues(e_w[i,:], time[i+1]-time[i]))
		# R_e[i+1,:,:] = R_e[i+1,:,:]/scipy.linalg.det(R_e[i+1,:,:])
		# print(R_e[i,:,:])

		w_d_dot[i,:] = (w_d[i+1,:]-w_d[i,:])/(time[i+1]-time[i])

	R[i,:,:] = np.dot(R_d[i,:,:], R_e[i,:,:])
	w[i,:] = e_w[i,:] + np.dot(R_e[i,:,:].transpose(), w_d[i,:])

	var = np.dot(hat(e_w[i,:]), np.dot(R_e[i,:,:].transpose(), w_d[i,:])) - np.dot(R_e[i,:,:].transpose(), w_d_dot[i,:])
	M_d[i,:] = -k_R*vee(e_R[i,:,:]) - k_w*e_w[i,:] + np.cross(w[i,:], np.dot(J,w[i,:])) - np.dot(J, var)

	Moment[i,:] =  M_e[i,:] + M_d[i,:]

	if i<N-1:
		Moment_dot[i+1,:] = (Moment[i+1]-Moment[i])/(time[i+1]-time[i])
		Moment_dot_dot[i+1,:] = (Moment_dot[i+1]-Moment_dot[i])/(time[i+1]-time[i])

	M_d_dot[i,:] = Moment_dot[i,:] - M_e_dot[i,:]

	if i<N-1:
		M_d_dot_dot[i+1,:] = (M_d_dot[i+1,:]-M_d_dot[i,:])/(time[i+1]-time[i])

	# print(Moment_dot[i,:])

	phi[i] = np.arcsin(R[i,2,1])
	th[i]  = np.arccos(R[i,2,2]/np.cos(phi[i]))
	psi[i] = np.arccos(R[i,1,1]/np.cos(phi[i]))

	T0 = m*g

	def del_int(t, x):
		return np.dot(np.array([[0.0, 1.0],[0.0, -2*x[1]*np.tan(x[0])]]), x) + np.array([0, Moment_dot_dot[i,2]/(2*l*T0)])

	def del_int_des(t, x):
		return np.dot(np.array([[0.0, 1.0],[0.0, -2*x[1]*np.tan(x[0])]]), x) + np.array([0, M_d_dot_dot[i,2]/(2*l*T0)])


	if i<N-1:
		x[i+1,:] = RK4(del_int, x[i], time[i], time[i+1], 0.001, 2)[0]
		y[i+1,:] = RK4(del_int_des, y[i], time[i], time[i+1], 0.001, 2)[0]
		delta[i+1] = x[i+1, 0]
		delta_dot[i+1] = x[i+1, 1]
		delta_des[i+1] = y[i+1,0]
		delta_des_dot[i+1] = y[i+1, 1]


	# delta[i] = np.arcsin(Moment[i,2]/(-2*l*T_m_norm))
	# delta_des[i] = np.arcsin(M_d[i,2]/(-2*l*T_m_norm))

	# print(delta[i])

	# print(M_d[i,2]/(-2*l*T_m_norm))


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
ax[1,1].plot(time, 180.0*(1/np.pi)*delta_des,'b')
ax[1,1].set_xlabel('Time')
ax[1,1].set_ylabel('delta')
ax[1,1].legend(['delta', 'des delta'], fontsize='xx-small')

plt.show()