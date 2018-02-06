from force_moments import *
from scipy import optimize as opt

class Trim(MAVForces):
    def __init(self):
        MAVForces.__init__(self)

        self.deltas = [0,0,0,0]
        self.update_init_conds = [0,0,0,0,0,0,0,0,0,0,0,0]

    def call_opt(self,star,x,trim,wind):
        self.Va_trim,self.gamma_trim,self.R_trim = trim
        minimum = opt.fmin(self.trim_opt,star,args=(x,trim,wind))

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

        self.a_theta1 = -0.5 * self.rho * self.Va**2 * self.c * self.S * self.C_mq * self.c / (self.Jy * 2 * self.Va_trim)
        self.a_theta2 = -0.5 * self.rho * self.Va**2 * self.c * self.S * self.C_malpha / self.Jy
        self.a_theta3 = 0.5 * self.rho * self.Va**2 * self.c * self.S * self.C_mdeltae / self.Jy

        self.a_v1 = self.rho * self.Va_trim * self.S / self.mass * (self.C_D0 + self.C_Dalpha * self.alpha_trim + self.C_Ddeltae * self.deltas[0]) + self.rho * self.S_prop / self.mass * self.C_prop * self.Va_trim
        self.a_v2 = self.rho * self.S_prop / self.mass * self.C_prop * self.k_motor**2 * self.deltas[3]
        self.a_v3 = self.g #* np.cos(th - self.chi)
        
        # Transfer Functions
        self.T_phi_deltaa = sigs.TransferFunction([self.a_phi2],[1,self.a_phi1,0])
        self.T_chi_phi = sigs.TransferFunction([self.g / self.Va_trim],[1,0])
        self.T_th_deltae = sigs.TransferFunction([self.a_theta3],[1,self.a_theta1,self.a_theta2])
        self.T_h_theta = sigs.TransferFunction([self.Va_trim],[1,0])
        self.T_h_Va = sigs.TransferFunction([th],[1,0])
        self.T_Va_deltat = sigs.TransferFunction([self.a_v2],[1,self.a_v1])
        self.T_Va_theta = sigs.TransferFunction([-self.a_v3],[1,self.a_v1])
        self.T_v_deltar = sigs.TransferFunction([self.Va_trim * self.a_beta2],[1,self.a_beta1])

        # print self.T_phi_deltaa,'\n'
        # print self.T_chi_phi,'\n'
        # print self.T_th_deltae,'\n'
        # print self.T_h_theta,'\n'
        # print self.T_h_Va,'\n'
        # print self.T_Va_deltat,'\n'
        # print self.T_Va_theta,'\n'
        # print self.T_v_deltar,'\n'

    def get_linearized_ss(self,x,wind):
        # Relabel Inputs
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        # fx,fy,fz,l,m,n,_,_,_,_,_,_ = self.force_calc(x,self.deltas,wind)
        # print m

        # Lateral SS Model Coefficients
        Yv = self.rho * self.S * self.b * v * (self.C_Yp * p + self.C_Yr * r) / (4 * self.mass * self.Va_trim) + self.rho * self.S * v / self.mass * (self.C_Y0 + self.C_Ybeta * self.beta + self.C_Ydeltaa * self.deltas[1] + self.C_Ydeltar * self.deltas[2]) + self.rho * self.S * self.C_Ybeta / (2 * self.mass) * np.sqrt(u**2 + w**2)
        Yp = w + self.rho * self.Va_trim * self.S * self.b / (4 * self.mass) * self.C_Yp
        Yr = -u + self.rho * self.Va_trim * self.S * self.b / (4 * self.mass) * self.C_Yr
        Ydeltaa = self.rho * self.Va_trim**2 * self.S / (2 * self.mass) * self.C_Ydeltaa
        Ydeltar = self.rho * self.Va_trim**2 * self.S / (2 * self.mass) * self.C_Ydeltar
        Lv = self.rho * self.S * self.b**2 * v / (4 * self.Va_trim) * (self.C_pp * p + self.C_pr * r) + self.rho * self.S * self.b * v * (self.C_p0 + self.C_pbeta * self.beta + self.C_pdeltaa * self.deltas[1] + self.C_pdeltar * self.deltas[2]) + self.rho * self.S * self.b * self.C_pbeta / 2 * np.sqrt(u**2 + w**2)
        Lp = self.g1 * q + self.rho * self.Va_trim * self.S * self.b**2 * self.C_pp / 4
        Lr = -self.g2 * q + self.rho * self.Va_trim * self.S * self.b**2 * self.C_pr / 4
        Ldeltaa = self.rho * self.Va_trim**2 * self.S * self.b * self.C_pdeltaa * 0.5
        Ldeltar = self.rho * self.Va_trim**2 * self.S * self.b * self.C_pdeltar * 0.5
        Nv = self.rho * self.S * self.b**2 * v / (4 * self.Va_trim) * (self.C_rp * p + self.C_rr * r) + self.rho * self.S * self.b * v * (self.C_r0 + self.C_rbeta * self.beta + self.C_rdeltaa * self.deltas[1] + self.C_rdeltar * self.deltas[2]) + self.rho * self.S * self.b * self.C_rbeta / 2 * np.sqrt(u**2 + w**2)
        Np = self.g7 * q + self.rho * self.Va_trim * self.S * self.b**2 * self.C_rp / 4
        Nr = -self.g1 * q + self.rho * self.Va_trim * self.S * self.b**2 * self.C_rr / 4
        Ndeltaa = self.rho * self.Va_trim**2 * self.S * self.b * self.C_rdeltaa * 0.5
        Ndeltar = self.rho * self.Va_trim**2 * self.S * self.b * self.C_rdeltar * 0.5

        # New coefficients needed
        self.C_X0 = -self.C_D0
        self.C_Xalpha = -self.C_Dalpha + self.C_L0
        self.C_Xq = self.calc_Cx_q(self.alpha_trim)
        self.C_Xdeltae = self.calc_Cx_deltae(self.alpha_trim)
        self.C_Z0 = -self.C_L0
        self.C_Zalpha = -self.C_D0 + self.C_Lalpha
        self.C_Zq = self.calc_Cz_q(self.alpha_trim)
        self.C_Zdeltae = self.calc_Cz_deltae(self.alpha_trim)
        
        # Longitudinal SS Model Coefficients
        Xu = u * self.rho * self.S / self.mass * (self.C_X0 + self.C_Xalpha * self.alpha_trim + self.C_Xdeltae * self.deltas[0]) - self.rho * self.S * w * self.C_Xalpha / (self.mass * 2) + self.rho * self.S * self.c * self.C_Xq * u * q / (4 * self.mass * self.Va_trim) - self.rho * self.S_prop * self.C_prop * u / self.mass
        Xw = -q + w * self.rho * self.S / self.mass * (self.C_X0 + self.C_Xalpha * self.alpha_trim + self.C_Xdeltae * self.deltas[0]) + self.rho * self.S * self.c * self.C_Xq * w * q / (4 * self.Va_trim * self.mass) + self.rho * self.S * self.C_Xalpha * u / (self.mass * 2) - self.rho * self.S_prop * self.C_prop * w / self.mass
        Xq = -w + self.rho * self.Va_trim * self.S * self.C_Xq * self.c / (4 * self.mass)
        Xdeltae = self.rho * self.Va_trim**2 * self.S * self.C_Xdeltae * self.c / (2 * self.mass)
        Xdeltat = self.rho * self.S_prop * self.C_prop * self.k_motor**2 * self.deltas[3] / self.mass
        Zu = q + u * self.rho * self.S / self.mass * (self.C_Z0 + self.C_Zalpha * self.alpha_trim + self.C_Zdeltae * self.deltas[0]) - self.rho * self.S * self.C_Zalpha * w / (2 * self.mass) + u * self.rho * self.S * self.C_Zq * self.c * q / (4 * self.mass * self.Va_trim)
        Zw = w * self.rho * self.S / self.mass * (self.C_Z0 + self.C_Zalpha * self.alpha_trim + self.C_Zdeltae * self.deltas[0]) + self.rho * self.S * self.C_Zalpha * u / (2 * self.mass) + self.rho * w * self.S * self.c * self.C_Zq * q / (4 * self.mass * self.Va_trim)
        Zq = u + self.rho * self.Va_trim * self.S * self.C_Zq * self.c / (4 * self.mass)
        Zdeltae = self.rho * self.Va_trim**2 * self.S * self.C_Zdeltae / (2 * self.mass)
        Mu = u * self.rho * self.S * self.c / self.Jy * (self.C_m0 + self.C_malpha * self.alpha_trim + self.C_mdeltae * self.deltas[0]) - self.rho * self.S * self.c * self.C_malpha * w / (2 * self.Jy) + self.rho * self.S * self.c**2 * self.C_mq * q * u / (4 * self.Jy * self.Va_trim)
        Mw = w * self.rho * self.S * self.c / self.Jy * (self.C_m0 + self.C_malpha * self.alpha_trim + self.C_mdeltae * self.deltas[0]) + self.rho * self.S * self.c * self.C_malpha * u / (2 * self.Jy) + self.rho * self.S * self.c**2 * self.C_mq * q * w / (4 * self.Jy * self.Va_trim)
        Mq = self.rho * self.Va_trim * self.S * self.c**2 * self.C_mq / (4 * self.Jy)
        Mdeltae = self.rho * self.Va_trim**2 * self.S * self.c * self.C_mdeltae / (2 * self.Jy)

        # Lateral State Space Matrices
        self.A_lat = [[Yv, Yp, Yr, self.g * np.cos(th) * np.cos(phi), 0],
                      [Lv, Lp, Lr, 0, 0],
                      [Nv, Np, Nr, 0, 0],
                      [0, 1, np.cos(phi) * np.tan(th), q * np.cos(phi) * np.tan(th) - r * np.sin(phi) * np.tan(th), 0],
                      [0, 0, np.cos(phi) * 1 / np.cos(th), p * np.cos(phi) * 1 / np.cos(th) - r * np.sin(phi) * 1 / np.cos(th), 0]]
        
        self.B_lat = [[Ydeltaa, Ydeltar],
                      [Ldeltaa, Ldeltar],
                      [Ndeltaa, Ndeltar],
                      [0, 0],
                      [0, 0]]

        # Longitudinal State Space matrices
        self.A_lon = [[Xu, Xw, Xq, -self.g * np.cos(th), 0],
                      [Zu, Zw, Zq, -self.g * np.sin(th), 0],
                      [Mu, Mw, Mq, 0, 0],
                      [0, 0, 1, 0, 0],
                      [np.sin(th), -np.cos(th), 0, u * np.cos(th) + w * np.sin(th), 0]]

        self.B_lon = [[Xdeltae, Xdeltat],
                      [Zdeltae, 0],
                      [Mdeltae, 0],
                      [0, 0],
                      [0, 0]]
                      
        
        
