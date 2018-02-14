from force_moments import *
from scipy import optimize as opt
from copy import deepcopy

class Trim(MAVForces):
    def __init(self):
        MAVForces.__init__(self)

    def variables(self):
        self.trim = [20.0,0.0,5e60]
        self.deltas = [0,0,0,0]
        self.update_init_conds = [0,0,0,0,0,0,0,0,0,0,0,0]
        self.star = [0.0,0.0,0.0]

    def call_opt(self,x,wind):

        self.Va_trim,self.gamma_trim,self.R_trim = self.trim
        minimum = opt.fmin(self.trim_opt,self.star,args=(x,self.trim,wind))
        self.delta_trim = deepcopy(self.deltas)
        return self.update_init_conds

    def trim_opt(self,star,x,trim,wind):
        ''' 
        Variables: star --> alpha, beta, phi
        States: x --> pn,pe,pd,u,v,w,phi,th,psi,p,q,r
        Desired Trim: --> Va,gamma,R
        Wind: --> wind
        '''
        # Relabel inputs
        alpha,beta,phi = star
        pn,pe,pd,u,v,w,phi_,th,psi,p,q,r = x
        Va,gamma,R = trim

        # Compute Trimmed States
        u = Va * np.cos(alpha) * np.cos(beta)
        v = Va * np.sin(beta)
        w = Va * np.sin(alpha) * np.cos(beta)
        th = alpha + gamma
        p = -Va / R * np.sin(th)
        q = Va / R * np.sin(phi) * np.cos(th)
        r = Va / R * np.cos(phi) * np.cos(th)

        # Compute Trimmed Input
        ele = (((self.Jx * (p**2 - r**2) + (self.Jx - self.Jz)*p*r)/(0.5 * self.rho * Va**2 * self.c * self.S)) - self.C_m0 - self.C_malpha * alpha - self.C_mq * self.c * q / (2 * Va)) /self.C_mdeltae
        
        thr = np.sqrt((2 * self.mass * (-r * v + q * w  + self.g * np.sin(th)) - self.rho * Va**2 * self.S * (self.calc_Cx(alpha) + self.calc_Cx_q(alpha) * self.c * q / (2 * Va) + self.calc_Cx_deltae(alpha) * ele)) / (self.rho * self.S_prop * self.C_prop * self.k_motor**2) + Va**2/(self.k_motor**2))

        C = [[self.C_pdeltaa, self.C_pdeltar],[self.C_rdeltaa,self.C_rdeltar]]
        C = np.linalg.inv(C)

        stuff = [(-self.g1 * p * q + self.g2 * q * r)/(0.5 * self.rho * Va**2 * self.S * self.b) - self.C_p0 - self.C_pbeta * beta - self.C_pp * self.b * p/(2 * Va) - self.C_pr * self.b * r/(2 * Va),
                 (-self.g7 * p * q + self.g1 * q * r)/(0.5 * self.rho * Va**2 * self.S * self.b) - self.C_r0 - self.C_rbeta * beta - self.C_rp * self.b * p/(2 * Va) - self.C_rr * self.b * r/(2 * Va)]

        delt = np.matmul(C,stuff)

        ail = delt[1]
        rud = delt[0]
        
        deltastar = [ele,ail,rud,thr]

        # Set to class variable
        self.deltas = deltastar
        self.update_init_conds = [pn,pe,pd,u,v,w,phi,th,psi,p,q,r]
        self.alpha_trim = alpha

        # Compute f(xstar,ustar)
        xdotstar = [0,0,Va * np.sin(gamma),0,0,0,0,0,Va / R * np.cos(gamma),0,0,0]
        xstar = [pn,pe,pd,u,v,w,phi,th,psi,p,q,r]
        fx,fy,fz,l,m,n,_,_,_,_,_,_ = self.force_calc(xstar,deltastar,wind)
        fxu = self.eom(xstar,0.1,fx,fy,fz,l,m,n)

        # Compute J
        J = np.linalg.norm(np.subtract(xdotstar,fxu))

        return J

    def get_transfer_functions(self,x):
        # Relabel Inputs
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        
        # Coefficients needed for transfer functions
        self.a_phi1 = -0.5 * self.rho * self.Va**2 * self.S * self.b * self.C_pp * self.b / (2 * self.Va)
        self.a_phi2 = 0.5 * self.rho * self.Va**2 * self.S * self.b * self.C_pdeltaa

        self.a_beta1 = -0.5 * self.rho * self.Va * self.S * self.C_Ybeta / self.mass
        self.a_beta2 = 0.5 * self.rho * self.Va * self.S * self.C_Ydeltar / self.mass

        self.a_th1 = -0.5 * self.rho * self.Va**2 * self.c * self.S * self.C_mq * self.c / (self.Jy * 2 * self.Va_trim)
        self.a_th2 = -0.5 * self.rho * self.Va**2 * self.c * self.S * self.C_malpha / self.Jy
        self.a_th3 = 0.5 * self.rho * self.Va**2 * self.c * self.S * self.C_mdeltae / self.Jy

        self.a_v1 = self.rho * self.Va_trim * self.S / self.mass * (self.C_D0 + self.C_Dalpha * self.alpha_trim + self.C_Ddeltae * self.deltas[0]) + self.rho * self.S_prop / self.mass * self.C_prop * self.Va_trim
        self.a_v2 = self.rho * self.S_prop / self.mass * self.C_prop * self.k_motor**2 * self.deltas[3]
        self.a_v3 = self.g #* np.cos(th - self.chi)
        
        # Transfer Functions
        self.T_phi_deltaa = sigs.TransferFunction([self.a_phi2],[1,self.a_phi1,0])
        self.T_chi_phi = sigs.TransferFunction([self.g / self.Va_trim],[1,0])
        self.T_th_deltae = sigs.TransferFunction([self.a_th3],[1,self.a_th1,self.a_th2])
        self.T_h_theta = sigs.TransferFunction([self.Va_trim],[1,0])
        self.T_h_Va = sigs.TransferFunction([th],[1,0])
        self.T_Va_deltat = sigs.TransferFunction([self.a_v2],[1,self.a_v1])
        self.T_Va_theta = sigs.TransferFunction([-self.a_v3],[1,self.a_v1])
        self.T_v_deltar = sigs.TransferFunction([self.Va_trim * self.a_beta2],[1,self.a_beta1])

    def get_numerical_ss(self,x):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        fx,fy,fz,l,m,n,_,_,_,_,_,_ = self.force_calc(x,self.deltas,[0,0,0,0,0,0])
        xdot_nom = self.eom(x,0.01,fx,fy,fz,l,m,n)
        xdot = np.zeros((12,12))
        step = 1e-8
        for i in xrange(12):
            for j in xrange(12):
                x0 = [pn,pe,pd,u,v,w,phi,th,psi,p,q,r]
                x0[j] += step
                fx,fy,fz,l,m,n,_,_,_,_,_,_ = self.force_calc(x0,self.deltas,[0,0,0,0,0,0])
                xdotij = self.eom(x0,0.01,fx,fy,fz,l,m,n)
                xdot[i,j] = (-xdot_nom[i] + xdotij[i]) / step
        self.A_lon = [[xdot[3,3],xdot[3,5],xdot[3,10],xdot[3,7],-xdot[3,2]],
                      [xdot[5,3],xdot[5,5],xdot[5,10],xdot[5,7],-xdot[5,2]],
                      [xdot[10,3],xdot[10,5],xdot[10,10],xdot[10,7],-xdot[10,2]],
                      [xdot[7,3],xdot[7,5],xdot[7,10],xdot[7,7],-xdot[7,2]],
                      [-xdot[2,3],-xdot[2,5],-xdot[2,10],-xdot[2,7],-xdot[2,2]]]

        self.A_lat = [[xdot[4,4],xdot[4,9],xdot[4,11],xdot[4,6],xdot[4,8]],
                      [xdot[9,4],xdot[9,9],xdot[9,11],xdot[9,6],xdot[9,8]],
                      [xdot[11,4],xdot[11,9],xdot[11,11],xdot[11,6],xdot[11,8]],
                      [xdot[6,4],xdot[6,9],xdot[6,11],xdot[6,6],xdot[6,8]],
                      [xdot[8,4],xdot[8,9],xdot[8,11],xdot[8,6],xdot[8,8]]]

        self.eig_lon,v = np.linalg.eig(self.A_lon)
        self.eig_lat,v = np.linalg.eig(self.A_lat)
        xdot = np.zeros((12,4))
        for i in xrange(12):
            for j in xrange(4):
                delta0 = deepcopy(self.deltas)
                delta0[j] += step
                fx,fy,fz,l,m,n,_,_,_,_,_,_ = self.force_calc(x,delta0,[0,0,0,0,0,0])
                xdotij = self.eom(x,0.01,fx,fy,fz,l,m,n)
                xdot[i,j] = (-xdot_nom[i] + xdotij[i]) / step
        self.B_lon = [[xdot[3,0],xdot[3,3]],
                      [xdot[5,0],xdot[5,3]],
                      [xdot[10,0],xdot[10,3]],
                      [xdot[7,0],xdot[7,3]],
                      [-xdot[2,0],-xdot[2,3]]]
        self.B_lat = [[xdot[4,1],xdot[4,2]],
                      [xdot[9,1],xdot[9,2]],
                      [xdot[11,1],xdot[11,2]],
                      [xdot[6,1],xdot[6,2]],
                      [xdot[8,1],xdot[8,2]]]
        
        
