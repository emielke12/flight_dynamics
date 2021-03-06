from airframe import AirCraftDrawing
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

class MAVEOM(AirCraftDrawing):
    def __init__(self):
        # Animation Class Inheritance
        AirCraftDrawing.__init__(self)

        # Time Step
        self.dt = 0.02
        self.t_sim = 0.0

        # Plane Parameters
        self.mass = 13.5
        self.g = 9.81
        self.Jx = 0.8244
        self.Jy = 1.135
        self.Jz = 1.759
        self.Jxz = 0.1204
        self.S = 0.55
        self.b = 2.8956
        self.c = 0.18994
        self.S_prop = 0.2027
        self.rho = 1.2682
        self.k_motor = 80.0
        self.k_Tp = 0.0
        self.k_Omega = 0.0
        self.e = 0.9
        self.C_L0 = 0.28
        self.C_D0 = 0.03
        self.C_m0 = -0.02338
        self.C_Lalpha = 3.45
        self.C_Dalpha = 0.3
        self.C_malpha = -0.38
        self.C_Lq = 0.0
        self.C_Dq = 0.0
        self.C_mq = -3.6
        self.C_Ldeltae = -0.36
        self.C_Ddeltae = 0.0
        self.C_mdeltae = -0.5
        self.C_prop = 1.0
        self.M = 50.0
        self.alpha_0 = 0.4712
        self.epsilon = 0.1592
        self.C_Dp = 0.0437
        self.C_Y0 = 0.0
        self.C_l0 = 0.0
        self.C_n0 = 0.0
        self.C_Ybeta = -0.98
        self.C_lbeta = -0.12
        self.C_nbeta = 0.25
        self.C_Yp = 0.0
        self.C_lp = -0.26
        self.C_np = 0.022
        self.C_Yr = 0.0
        self.C_lr = 0.14
        self.C_nr = -0.35
        self.C_Ydeltaa = 0.0
        self.C_ldeltaa = 0.08
        self.C_ndeltaa = 0.06
        self.C_Ydeltar = -0.17
        self.C_ldeltar = 0.105
        self.C_ndeltar = -0.032

        g = self.Jx * self.Jz - self.Jxz**2
        self.g1 = self.Jxz * (self.Jx - self.Jy + self.Jz) / g
        self.g2 = (self.Jz * (self.Jz - self.Jy) + self.Jxz**2) / g
        self.g3 = self.Jz / g
        self.g4 = self.Jxz / g
        self.g5 = (self.Jz - self.Jx) / self.Jy
        self.g6 = self.Jxz / self.Jy
        self.g7 = ((self.Jx - self.Jy) * self.Jx + self.Jxz**2)/g
        self.g8 = self.Jx / g

        self.C_p0 =  self.g3 * self.C_l0 + self.g4 * self.C_n0
        self.C_pbeta = self.g3 * self.C_lbeta + self.g4 * self.C_nbeta
        self.C_pp = self.g3 * self.C_lp + self.g4 * self.C_np
        self.C_pr = self.g3 * self.C_lr + self.g4 * self.C_nr
        self.C_pdeltaa = self.g3 * self.C_ldeltaa + self.g4 * self.C_ndeltaa
        self.C_pdeltar = self.g3 * self.C_ldeltar + self.g4 * self.C_ndeltar
        self.C_r0 = self.g4 * self.C_l0 + self.g8 * self.C_n0
        self.C_rbeta = self.g4 * self.C_lbeta + self.g8 * self.C_nbeta
        self.C_rp = self.g4 * self.C_lp + self.g8 * self.C_np
        self.C_rr = self.g4 * self.C_lr + self.g8 * self.C_nr
        self.C_rdeltaa = self.g4 * self.C_ldeltaa + self.g8 * self.C_ndeltaa
        self.C_rdeltar = self.g4 * self.C_ldeltar + self.g8 * self.C_ndeltar


    def eom(self,x,t,fx,fy,fz,l,m,n):
        # States
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        
        # Eq 3.14 Translational Kinematics
        R = [[np.cos(th) * np.cos(psi), np.sin(phi) * np.sin(th) * np.cos(psi) - np.cos(phi) * np.sin(psi), np.cos(phi) * np.sin(th) * np.cos(psi) + np.sin(phi) * np.sin(psi)],
             [np.cos(th) * np.sin(psi), np.sin(phi) * np.sin(th) * np.sin(psi) + np.cos(phi) * np.cos(psi), np.cos(phi) * np.sin(th) * np.sin(psi) - np.sin(phi) * np.cos(psi)],
             [-np.sin(th), np.sin(phi) * np.cos(th), np.cos(phi) * np.cos(th)]]
        dpNEDdt = np.matmul(R,[u,v,w])
        
        # Eq 3.15 Translational Dynamics
        duvwdt = [r * v - q * w + fx/self.mass,
                  p * w - r * u + fy/self.mass,
                  q * u - p * v + fz/self.mass]

        # Eq 3.16 Rotational Kinematics
        R = [[1, np.sin(phi)*np.tan(th), np.cos(phi)*np.tan(th)],
             [0, np.cos(phi), -np.sin(phi)],
             [0, np.sin(phi)/np.cos(th), np.cos(phi)/np.cos(th)]]
        deuldt = np.matmul(R,[p,q,r])

        # Eq 3.17 Rotational Dynamics
        dpqrdt = [self.g1 * p * q - self.g2 * q * r + self.g3 * l + self.g4 * n,
                  self.g5 * p * r - self.g6 * (p**2 - r**2) + m / self.Jy,
                  self.g7 * p * q - self.g1 * q * r + self.g4 * l + self.g8 * n]

        # Total Equation of motion outputs
        xdot = [dpNEDdt[0], dpNEDdt[1], dpNEDdt[2], duvwdt[0], duvwdt[1], duvwdt[2], 
                deuldt[0], deuldt[1], deuldt[2], dpqrdt[0], dpqrdt[1], dpqrdt[2]]
        return xdot

