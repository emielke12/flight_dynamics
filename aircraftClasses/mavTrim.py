from force_moments import *
from scipy import optimize as opt

class Trim(MAVForces):
    def __init(self):
        MAVForces.__init__(self)

        self.deltas = [0,0,0,0]
        self.update_init_conds = [0,0,0,0,0,0,0,0,0,0,0,0]

    def call_opt(self,star,x,trim,wind):
        minimum = opt.fmin(self.trim_opt,star,args=(x,trim,wind))

    def trim_opt(self,star,x,trim,wind):
        ''' 
        Variables: star --> alpha, beta, phi
        States: x --> u,v,w,th,p,q,r
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

        ail = delt[0]
        rud = delt[1]
        
        deltastar = [ele,ail,rud,thr]

        # Set to class variable
        self.deltas = deltastar
        self.update_init_conds = [pn,pe,pd,u,v,w,phi,th,psi,p,q,r]

        # Compute f(xstar,ustar)
        xdotstar = [0,0,Va * np.sin(gamma),0,0,0,0,0,Va / R * np.cos(gamma),0,0,0]
        xstar = [pn,pe,pd,u,v,w,phi,th,psi,p,q,r]
        fx,fy,fz,l,m,n,n1,n2,n3,n4,n5,n6 = self.force_calc(xstar,deltastar,wind)
        fxu = self.eom(xstar,0.1,fx,fy,fz,l,m,n)

        # Compute J
        J = np.linalg.norm(np.subtract(xdotstar,fxu))

        return J
