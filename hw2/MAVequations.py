from airframe import AirCraftDrawing
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class MAVEOM(AirCraftDrawing):
    def __init__(self):
        # Animation Class Inheritance
        AirCraftDrawing.__init__(self)

        # Plane Parameters
        self.M = 13.5
        self.g = 9.81
        self.Jx = 0.8244
        self.Jy = 1.135
        self.Jz = 1.759
        self.Jxz = 0.1204

    def eom(self,x,t,fx,fy,fz,l,m,n):
        # States
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        
        # Eq 3.14 Translational Kinematics
        R = np.transpose([[np.cos(th)*np.cos(psi), np.cos(th)*np.sin(psi), -np.sin(th)],
             [np.sin(phi)*np.sin(th)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(th)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(th)],
             [np.cos(phi)*np.sin(th)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(th)]])
        dpNEDdt = np.matmul(R,[u,v,w])
        
        # Eq 3.15 Translational Dynamics
        R = [r * v - q * w,
             p * w - r * u,
             q * u - p * v]
        duvwdt = np.multiply(1.0/self.M,[fx,fy,fz]) + R

        # Eq 3.16 Rotational Kinematics
        R = [[1, np.sin(phi)*np.tan(th), np.cos(phi)*np.tan(th)],
             [0, np.cos(phi), -np.sin(phi)],
             [0, np.sin(phi)/np.cos(th), np.cos(phi)/np.cos(th)]]
        deuldt = np.matmul(R,[p,q,r])

        # Eq 3.17 Rotational Dynamics
        g = self.Jx*self.Jz - self.Jxz**2
        g1 = self.Jxz*(self.Jx - self.Jy + self.Jz)
        g2 = self.Jz * (self.Jz - self.Jy) + self.Jxz**2
        g3 = self.Jz / g
        g4 = self.Jxz / g
        g5 = (self.Jz - self.Jx) / self.Jy
        g6 = self.Jxz / self.Jy
        g7 = ((self.Jx - self.Jy) * self.Jx + self.Jxz**2)/g
        g8 = self.Jx / g

        R1 = [g1 * p * q - g2 * q * r,
              g5 * p * r - g6 * (p**2 - r**2),
              g7 * p * q - g1 * q * r]
        R2 = [g3 * l + g4 * n,
              m / self.Jy,
              g4 * l + g8 * n]
        dpqrdt = np.add(R1,R2).tolist()

        xdot = [dpNEDdt[0], dpNEDdt[1], dpNEDdt[2], duvwdt[0], duvwdt[1], duvwdt[2], 
                deuldt[0], deuldt[1], deuldt[2], dpqrdt[0], dpqrdt[1], dpqrdt[2]]
        return xdot

if __name__ == "__main__":
    # aircraft = AirCraftDrawing()
    meq = MAVEOM()

    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = 0
    u0 = 0
    v0 = 0
    w0 = 0
    phi0 = 0
    th0 = 0
    psi0 = 0
    p0 = 0
    q0 = 0
    r0 = 0

    x0 = [pn0,pe0,pd0,u0,v0,w0,phi0,th0,psi0,p0,q0,r0]

    # Time Steps
    step = 301

    # Inputs
    fx = np.zeros(step)
    fy = np.zeros(step)
    fz = np.zeros(step)
    l = np.zeros(step)
    m = np.zeros(step)
    n = np.zeros(step)
    # n[50:70] = 10.0
    # n[70:110] = -10.0
    # n[110:130] = 10.0

    # Time Vector
    t = np.linspace(0,10,step)

    # Solve ODE
    sols = np.zeros((len(t),len(x0)))
    for i in xrange(step):
        sol = odeint(meq.eom,x0,[t[i-1],t[i]],args=(fx[i],fy[i],fz[i],l[i],m[i],n[i]))
        for j in xrange(len(x0)):
            sols[i,j] = sol[-1][j]
        x0 = sol[-1:][0]        


    pos = [sols[:,0],sols[:,1],sols[:,2]]
    eul = [sols[:,6],sols[:,7],sols[:,8]]
    meq.input_goal(eul,pos)
