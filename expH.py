
# coding: utf-8

# In[4]:


import numpy as np
import scipy.sparse.linalg
import fermion as f

dim = 4
spin = 2
D = spin*dim
h = f.fermion(dim, spin)
#hopping matrix(2 by2 nearest neighbor)
t = -1.
t_mat = np.array([[0., t, 0., t],
                  [t, 0., t, 0.],
                  [0., t, 0., t],
                  [t, 0., t, 0.]])

zero = np.zeros((dim,dim))
t = np.insert(t_mat,range(1,1+dim), zero, axis=1)
s = np.array([[t[j, i-1] for i in range(dim*2)]for j in range(dim)])
t_mat = np.insert(t, range(1, 1+dim), s, axis=0)

#interaction (on site)
U_mat = np.zeros((D**2, D**2))
U=3.
l = np.zeros((D,D))
for i in range(dim):
    l[2*i, 2*i+1] = 1
l = l.T + l
l = l.reshape(D**2)
for i in np.where(l==1)[0]:
    U_mat[i,i] = -U/2

h1 = h.h1(t_mat)
h2 = h.h2(U_mat)
ham = h1 + h2

def Ham(n):
    l = h.n_sector(n)
    return f.part(ham, l, l).tocsc()
def gibbs_state(beta, n):
    vec = scipy.sparse.linalg.expm_multiply(- beta*Ham(n) ,np.random.rand(Ham(n).shape[0])+1j*np.random.rand(Ham(n).shape[0]))
    return vec

v = gibbs_state(50, 4)
np.dot(np.conjugate(v), Ham(4)@v)/np.dot(np.conjugate(v), v)


# In[33]:


def energy(beta, n):
    v = gibbs_state(beta, n)
    return np.dot(np.conjugate(v), Ham(n)@v)/np.dot(np.conjugate(v),v)

x = np.array([energy(50, 4) for i in range(10)])
np.average(x)

