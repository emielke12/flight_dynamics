from airframe import AirCraftDrawing
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class MAVEOM(AirCraftDrawing):
    def __init__(self):
        # Animation Class Inheritance
        AirCraftDrawing.__init__(self)

        # Plane Parameters
        self.M = 2
        self.g = 9.81
        self.Jx = 0.01
        self.Jy = 0.01
        self.Jz = 0.01
        self.Jxz = 0.01

    def equations_of_motion(self,x,t,uu):
        # States
        pn = x[0]
        pe = x[1]
        pd = x[2]
        u = x[3]
        v = x[4]
        w = x[5]
        phi = x[6]
        th = x[7]
        psi = x[8]
        p = x[9]
        q = x[10]
        r = x[11]

        # Inputs
        fx = uu[0]
        fy = uu[1]
        fz = uu[2]
        l = uu[3]
        m = uu[4]
        n = uu[5]
        
        # Eq 3.14 Translational Kinematics
        R = np.transpose([[np.cos(th)*np.cos(psi), np.cos(th)*np.sin(psi), -np.sin(th)],
             [np.sin(phi)*np.sin(th)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(th)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(th)],
             [np.cos(phi)*np.sin(th)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(th)]])
        dpNEDdt = np.matmul(R,[u,v,w])
        
        # Eq 3.15 Translational Dynamics
        R = [[r * v - q * w],
             [p * w - r * u],
             [q * u - p * v]]
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

        R1 = [[g1 * p * q - g2 * q * r],
              [g5 * p * r - g6 * (p**2 - r**2)],
              [g7 * p * q - g1 * q * r]]
        R2 = [[g3 * l + g4 * n],
              [m / self.Jy],
              [g4 * l + g8 * n]]
        dpqrdt = np.add(R1,R2)

        xdot = [dpNEDdt[0], dpNEDdt[1], dpNEDdt[2], duvwdt[0], duvwdt[1], duvwdt[2], 
                deuldt[0], deuldt[1], deuldt[2], dpqrdt[0], dpqrdt[1], dpqrdt[2]]
        return xdot


if __name__ == "__main__":
    aircraft = AirCraftDrawing()
    meq = MAVEOM(AirCraftDrawing)

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

    # Inputs
    fx = 1.0
    fy = 0.0
    fz = 0.0
    l = 0.0
    m = 0.0
    n = 0.0

    uu = [fx,fy,fz,l,m,n]

    # Time Vector
    t = np.linspace(0,10,101)

    # Solve ODE
    sol = odeint(equations_of_motion,x0,t,args=(uu))

#     # Plot
#     plt.plot(t,sol[:,0])
#     plt.show()
    pos = [sol[:,0],sol[:,1],sol[:,2]]
    eul = [sol[:,6],sol[:,7],sol[:,8]]
    meq.input_goal(eul,pos)
