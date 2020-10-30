import numpy as np
import scipy.linalg


def hat(v):
	A = np.zeros((3,3))
	A[0,1] = -v[2]
	A[0,2] = v[1]
	A[1,2] = -v[0]
	A[1,0] = v[2]
	A[2,0] = -v[1]
	A[2,1] = v[0]
	return A

def rodrigues(t, omega):
	if omega.all() == 0:
		K = np.zeros((3,3))
	else:
		K = hat(omega/scipy.linalg.norm(omega))
	th = scipy.linalg.norm(omega)*t
	return np.eye(3) + 	np.sin(th)*K + (1-np.cos(th))*np.dot(K,K)
print(rodrigues(1.0, np.array([1000000000.0, 10000000000000.0,10000.0])))
