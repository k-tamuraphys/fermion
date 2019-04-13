
# coding: utf-8

# In[90]:


import numpy as np
from scipy import sparse
#%debug
#Pauli matrices
p1 = sparse.csc_matrix(np.array([[0., 1.],
                   [1., 0.]]))
p2 = sparse.csc_matrix(np.array([[0., -1.j],
                   [1.j, 0.]]))
p3 = sparse.csc_matrix(np.array([[1., 0.],
                   [0., -1.]]))
raising_op = 1/2 *(p1 + 1.j*p2)
lowering_op = 1/2 *(p1 - 1.j*p2)


def vm_dot(v, m):#m = [m0, m1, ... m_n-1], v = (v0, ... , vn-1)
        dim = len(m)
        x = np.array([v[i]*m[i] for i in range(dim)])
        return sum(x)
    
def linear_tf(t, m):#t: n by n matrix, m: n dimensional vector of d by d matrix
    dim = len(m)
    x = np.array([vm_dot(t[i], m) for i in range(dim)])
    return x

def mm_dot(m1, m2):
    dim = len(m1)
    x = np.array([m1[i]@m2[i] for i in range(dim)])
    return sum(x)
    
def quad(m1, t, m2):
    m = linear_tf(t,m2)
    return mm_dot(m1, m)

def ndiag(d):
    a = np.array([0])
    for i in range(d):
        b = np.ones(2**i) + a
        a = np.append(b, a)
    return a

def part(a, l1, l2):
        L1 = np.array([l1]).T
        L2 = [l2]
        for k in range(len(l1)-1):
            L2 = np.insert(L2, 0, l2, axis=0)
        for k in range(len(l2)-1):
            L1 = np.insert(L1, 0, l1, axis=1)
        return a[L1, L2]
class ops:
    
    def vm_dot(self, v, m):#m = [m0, m1, ... m_n-1], v = (v0, ... , vn-1)
        dim = len(m)
        x = np.array([v[i]*m[i] for i in range(dim)])
        return sum(x)
    
    def linear_tf(self, t, m):#t: n by n matrix, m: n dimensional vector of d by d matrix
        dim = len(m)
        x = np.array([self.vm_dot(t[i], m) for i in range(dim)])
        return x

    def mm_dot(self, m1, m2):
        dim = len(m1)
        x = np.array([m1[i]@m2[i] for i in range(dim)])
        return sum(x)
    
    def quad(self, m1, t, m2):
        m = self.linear_tf(t,m2)
        return self.mm_dot(m1, m)
    
    def part(self, a, l1, l2):
        L1 = np.array([l1]).T
        L2 = [l2]
        for k in range(len(l1)-1):
            L2 = np.insert(L2, 0, l2, axis=0)
        for k in range(len(l2)-1):
            L1 = np.insert(L1, 0, l1, axis=1)
        return a[L1, L2]



class fermion:
    def __init__(self, dim, spin=2):
        self.dim = dim
        self.spin = spin
        
    def c(self, i, alpha):
        s = self.spin
        d = s*i+alpha
        c = 1.
        for j in range(d):
            c = sparse.kron(c, -p3, format="csc")
        
        c = sparse.kron(c, lowering_op, format="csc")
    
        for j in range(s*self.dim - 1 - d):
            c = sparse.kron(c, sparse.identity(2, format="csc"), format="csc")
        return c.tocsc()
    
    def cdag(self, i, alpha):
        return np.conjugate(self.c(i,alpha).T).tocsc()
    
    def c_vec(self):#(c(0,alpha=0), c(0, alpha=1), c(1, alpha=0), ... ,c(dim-1, 1))という順序で並んだfermion op.からなるベクトル
        return np.array([self.c(i,alpha) for i in range(self.dim) for alpha in range(self.spin)])
    
    def c_vec2(self):
        return np.array([self.c(i,alpha)@self.c(j,beta) for i in range(self.dim) for alpha in range(self.spin) for j in range(self.dim) for beta in range(self.spin)])
    
    def cdag_vec(self):
        return np.array([self.cdag(i,alpha) for i in range(self.dim) for alpha in range(self.spin)])
    
    def cdag_vec2(self):
        return np.array([self.cdag(i,alpha)@self.cdag(j,beta) for i in range(self.dim) for alpha in range(self.spin) for j in range(self.dim) for beta in range(self.spin)])
    
    def n_vec(self):
        return np.array([self.cdag(i,alpha)@self.c(i,alpha) for i in range(self.dim) for alpha in range(self.spin)])
    
    def n_tot(self):
        return sum(self.nvec())
    
    def n_sector(self, n):
        Ndiag = ndiag(self.spin*self.dim)
        l = np.where(Ndiag == n)[0]
        return l
    
    def h1(self, t_mat):
        return quad(self.cdag_vec(), t_mat, self.c_vec())
    
    def h2(self, U_mat):
        return quad(self.cdag_vec2(), U_mat, self.c_vec2())

