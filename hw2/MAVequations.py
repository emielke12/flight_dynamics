from airframe import AirCraftDrawing
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class MAVEOM(AirCraftDrawing):
    def __init__(self):
        M = 2

    def trans_kin(self,u,t,psi,th,phi):
        ''' Eq. 3.14 '''
        R = np.transpose([[np.cos(th)*np.cos(psi), np.cos(th)*np.sin(psi), -np.sin(th)],
             [np.sin(phi)*np.sin(th)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(th)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(th)],
             [np.cos(phi)*np.sin(th)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(th)]])
        dpedt = np.matmul(R,u)
        return dpedt

    def rot_kin(self,p,t,psi,th,phi):
        ''' Eq. 3.16 '''
        R = [[1, np.sin(phi)*np.tan(th), np.cos(phi)*np.tan(th)],
             [0, np.cos(phi), -np.sin(phi)],
             [0, np.sin(phi)/np.cos(th), np.cos(phi)/np.cos(th)]]
        deuldt = np.matmul(R,p)
        return deuldt

    def trans_dyn(self,f,t,p,q,r,u,v,w,M):
        ''' Eq. 3.15 '''
        R = [[r * v - q * w],
             [p * w - r * u],
             [q * u - p * v]]
        dudt = np.multiply(1.0/M,f) + R
        return dudt

    def rot_dyn(self,tau,t,p,q,r,Jx,Jy,Jz,Jxz):
        ''' Eq. 3.17 '''
        l,m,n = tau
        g = Jx*Jz - Jxz**2
        g1 = Jxz*(Jx - Jy + Jz)
        g2 = Jz * (Jz - Jy) + Jxz**2
        g3 = Jz / g
        g4 = Jxz / g
        g5 = (Jz - Jx) / Jy
        g6 = Jxz / Jy
        g7 = ((Jx - Jy) * Jx + Jxz**2)/g
        g8 = Jx / g

        R1 = [[g1 * p * q - g2 * q * r],
              [g5 * p * r - g6 * (p**2 - r**2)],
              [g7 * p * q - g1 * q * r]]
        R2 = [[g3 * l + g4 * n],
              [m / Jy],
              [g4 * l + g8 * n]]
        dpdt = np.add(R1,R2)
        return dpdt

# psi = 0.0
# th = 0.0
# phi = 0.0

# y0 = [1.0, 0.0, 0.0]

# t = np.linspace(0,10,101)

# sol = odeint(pend,y0,t,args=(psi,th,phi))

# plt.plot(t,sol[:,0])
# plt.show()
# B = MAVEOM()
# print B.fuse_l1
