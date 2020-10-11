import numpy as np
import scipy
import matplotlib.pyplot as plt
import RK4

################### CONSTANTS ############################
L1 = 0.8
L2 = 0.5
O1 = -0.666
O2 = 0.333
a1 = 0.04
a2 = 0.04
d  = 0.04
c  = 0.1
psi = 1
mu = 1
w_r = 0.25
a = 0.5
theta = np.pi/4.0

################### INITIALISE TIME #######################
end = 50
N = 10*end+1
t = np.linspace(0, end, N)

################### INITIALISE q ##########################
q= np.zeros((2,N))
q[:,0] = np.array([1.3, -1.3])

q_dot = np.zeros((2,N))
q_dot[:,0] = np.array([0.0, 0.0])

qD= np.zeros((2,N))
qD_dot = np.zeros((2,N))


################### Desired yD generation #################
r = np.zeros(2)

yD = np.zeros((2,N))
#yD[:,0] = np.array([0.010,0.016])
xD = np.zeros((2,N))
xD_dot = np.zeros((2,N))
xD_ddot = np.zeros((2,N))

################# Torque values ############################
tau_arr = np.zeros((2, N))
tau_star_arr = np.zeros((2, N))

################### Parameters #############################
p1 = 3.31
p2 = 0.116
p3 = 0.16
p = np.array([p1, p2, p3])

def yD_dot(t, yD):
	r1 = a1*np.sin(w_r*t) + c + d*np.sin(1.5*w_r*t)
	r2 = a2*np.sin(w_r*t + psi) + c + d*np.sin(1.5*w_r*t + psi)
	r = np.array([r1, r2])
	return -mu*yD + r

def r(t):
	return np.array([a1*np.sin(w_r*t) + c + d*np.sin(1.5*w_r*t), a2*np.sin(w_r*t + psi) + c + d*np.sin(1.5*w_r*t + psi)])

def r_dot(t):
	r1 = a1*w_r*np.cos(w_r*t) + 1.5*d*w_r*np.cos(1.5*w_r*t)
	r2 = a2*w_r*np.cos(w_r*t + psi) + 1.5*d*w_r*np.cos(1.5*w_r*t + psi)
	return np.array([r1, r2])

def omega(q):
	x1 = L1*np.cos(q[0]) + L2*np.cos(q[0]+q[1]) + O1
	x2 = L1*np.sin(q[0]) + L2*np.sin(q[0]+q[1]) + O2
	return np.array([x1, x2])

def J(q):
	x1 = -L1*np.sin(q[0]) - L2*np.sin(q[0]+q[1])
	x2 = -L2*np.sin(q[0]+q[1])
	x3 = L1*np.cos(q[0]) + L2*np.cos(q[0]+q[1])
	x4 = L2*np.cos(q[0]+q[1])
	return np.array([[x1, x2],[x3, x4]])

def J_inv(q):
	x1 = L2*np.cos(q[0]+q[1])
	x2 = L2*np.sin(q[0]+q[1])
	x3 = -L1*np.cos(q[0]) - L2*np.cos(q[0]+q[1])
	x4 = -L1*np.sin(q[0]) - L2*np.sin(q[0]+q[1])
	b = 1.0/(L1*L2*np.sin(q[1]))
	return b*np.array([[x1, x2],[x3, x4]])

def J_dot(q, q_dot):
	x1 = -L1*q_dot[0]*np.cos(q[0]) - L2*(q_dot[0]+q_dot[1])*np.cos(q[0]+q[1])
	x2 = -L2*(q_dot[0]+q_dot[1])*np.cos(q[0]+q[1])
	x3 = -L1*q_dot[0]*np.sin(q[0]) - L2*(q_dot[0]+q_dot[1])*np.sin(q[0]+q[1])
	x4 = -L2*(q_dot[0]+q_dot[1])*np.sin(q[0]+q[1])
	return np.array([[x1, x2],
					 [x3, x4]])

def M(q, p):
	x1 = p[0] + 2*p[2]*np.cos(q[1])
	x2 = p[1] + p[2]*np.cos(q[1])
	x3 = x2
	x4 = p[2]
	return np.array([[x1, x2],
					 [x3, x4]])

def V(q, q_dot, p):
	z = p[2]*np.sin(q[1])
	x1 = -q_dot[1]
	x2 = q_dot[0] + q_dot[1]
	x3 = q_dot[1]
	x4 = 0.0
	return z*np.array([[x1, x2],
					   [x3, x4]])


################### INITIALISE x and x_dot ################
x = np.zeros((2,N))
x[:,0] = omega(q[:,0])

x_dot = np.zeros((2,N))
x_dot[:,0] = np.dot(J(q[:,0]), q_dot[:,0])

z = np.zeros((4, N))


A = np.array([[np.cos(theta), np.sin(theta)],
			 [-np.sin(theta), np.cos(theta)]])


Kp = np.array([[1.0, 0.0],
			   [0.0, 1.0]])

# r = np.array([1,2,3])
# l = np.array([5,6,7])
# print(np.concatenate((r,l)))

for i in range(N-1):
	yD[:,i+1] = yD[:,i] + 0.1*yD_dot(t[i],yD[:,i])
	xD[:,i] = a*np.dot(A, yD[:,i])
	xD_dot[:,i] = a*np.dot(A, -mu*yD[:,i]+r(t[i]))
	xD_ddot[:,i] = a*np.dot(A, mu*mu*yD[:,i] - mu*r(t[i]) + r_dot(t[i]))
	M_star = np.dot(np.transpose(J_inv(q[:,i])), np.dot(M(q[:,i], p), J_inv(q[:,i])))
	V_star_st1 = np.dot(np.transpose(J_inv(q[:,i])), np.dot(V(q[:,i], q_dot[:,i], p), J_inv(q[:,i])))
	V_star = V_star_st1 - np.dot(M_star, np.dot(J_dot(q[:,i], q_dot[:,i]), J_inv(q[:,i])))

	e = x[:,i] - xD[:,i]
	e_dot = x_dot[:,i] - xD_dot[:,i]
	s = e + e_dot

	W1phi = np.dot(M_star, (xD_ddot[:,i] - e_dot)) + np.dot(V_star, (xD_dot[:,i] - e))
	tau_star = -np.dot(Kp, s) + W1phi
	tau_star_arr[:,i] = tau_star

	def int_dyn(t, z):
		A = np.zeros((4,4))
		A[0:2, 2:4] = np.eye(2)
		A[2:4, 0:2] = -(V_star + Kp)
		A[2:4, 2:4] = -(M_star + V_star + Kp)
		c1 = np.dot(M_star, xD_ddot[:,i]) + np.dot((M_star+V_star+Kp), xD_dot[:,i]) + np.dot((V_star+Kp), xD[:,i])
		d = np.zeros(4)
		d[2:4] = c1

		return np.dot(A, z) + d

	z[:,i] = np.concatenate((x[:,i], x_dot[:,i]))
	z[:,i+1] = z[:,i] + 0.1*int_dyn(t, z[:,i])

	x_dot[:, i+1] = z[:, i+1][2:4]
	x[:, i+1] = z[:, i+1][0:2]

	cosq2 = 0.5*((x[0,i+1]-O1)**2 + (x[1,i+1]-O2)**2 - (L1)**2 - (L2)**2)
	q[1, i+1] = np.arccos(cosq2)

	cosq2d = 0.5*((xD[0,i]-O1)**2 + (xD[1,i]-O2)**2 - (L1)**2 - (L2)**2)
	qD[1,i] = np.arccos(cosq2d)


	mag = np.sqrt(L1**2 + L2**2 + 2.0*L1*L2*cosq2)
	val = L1 + L2*cosq2
	th = np.arccos((L1 + L2*cosq2)/mag)

	q[0, i+1] = np.arcsin((x[1,i+1] - O2)/mag) - th

	magd = np.sqrt(L1**2 + L2**2 + 2.0*L1*L2*cosq2d)
	vald = L1 + L2*cosq2d
	thd = np.arccos((L1 + L2*cosq2d)/magd)

	qD[0, i] = np.arcsin((xD[1,i] - O2)/magd) - thd



fig, ax = plt.subplots(2,2)
ax[0,0].plot(t, 180.0*(1/np.pi)*q[0,:],'r')
ax[0,0].plot(t, 180.0*(1/np.pi)*qD[0,:],'b')
ax[0,0].set_xlabel('Time')
ax[0,0].set_ylabel('q1')
ax[0,0].legend(['q1', 'des pq1'], fontsize='xx-small')

ax[0,1].plot(t, 180.0*(1/np.pi)*q[1,:],'r')
ax[0,1].plot(t, 180.0*(1/np.pi)*qD[1,:],'b')
ax[0,1].set_xlabel('Time')
ax[0,1].set_ylabel('q2')
ax[0,1].legend(['q2', 'des q2'], fontsize='xx-small')

plt.show()




	



# plt.plot(t, 180.0*(1.0/np.pi)*q[0,:], 'b')
# plt.plot(t, 180.0*(1.0/np.pi)*q[1,:], 'r')
# plt.show()


	