import numpy as np

def power(A):
  x = A[0]
  for i in range(20):
    x = A @ x
    x /= np.linalg.norm(x)
  return x

def svd(A):
  nrow, ncol = A.shape
  if ncol > nrow: 
    A = A.T
    nrow, ncol = A.shape
  U = np.zeros((nrow, ncol))
  V = np.zeros((ncol, ncol))
  S = np.zeros(ncol)
  
  for i in range(ncol):
    v = power(A.T @ A)
    u = A @ v
    u = u/np.linalg.norm(u)
    s = u @ A @ v
    #print(u, v, s)
    A = A - s * u[:, None] * v
    U[:,i] = u
    V[:,i] = v
    S[i] = s
  return  (S, V, U)if(ncol > nrow) else (S, U, V)

# 
# A = np.random.rand(4000, 300)
# A = np.arange(1,5).astype(float).reshape(2,-1)
# a = svd(A)
# 
# 
# S = [[ 35.1826483319, 1.4769077000, 0.0000000000 ]]
# 
# 
# U = [[ 0.1013459963, -0.7679381414, 0.1848874205 ],
#  [ 0.2485687533, -0.4880712805, 0.2516837879 ],
#  [ 0.3957915103, -0.2082044196, 0.3110371519 ],
#  [ 0.5430142673, 0.0716624412, 0.4010830747 ],
#  [ 0.6902370243, 0.3515293021, 0.8030318611 ]]
# 
# 
# V = [[ 0.5192726084, 0.7507924423, 0.4146734714 ],
#  [ 0.5755207207, 0.0459263913, 0.7475281103 ],
#  [ 0.6317688329, -0.6589396597, 0.5188907750 ]]
