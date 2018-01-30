from MAVequations import *
import scipy.signal as sigs

class MAVForces(MAVEOM):
    def __init__(self):
        MAVEOM.__init__(self)
        self.Va = 0.0
        self.alpha = 0.0
        self.beta = 0.0
        
        self.low_altitude_low_turbulence = [50,200,50,1.06,0.7]
        self.low_altitude_mod_turbulence = [50,200,50,2.12,1.4]        
        self.med_altitude_low_turbulence = [600,533,533,1.5,1.5]        
        self.med_altitude_mod_turbulence = [600,533,533,3.0,3.0]
        self.gust = {'low_low':self.low_altitude_low_turbulence,
                     'low_mod':self.low_altitude_mod_turbulence,
                     'med_low':self.med_altitude_low_turbulence,
                     'med_mod':self.med_altitude_mod_turbulence}

    def force_calc(self,x, delta, wind):
        ''' Calclates forces and moments acting on airframe
        Inputs are states, control surface deltas, wind
        Outputs are Forces, torques, wind speed, alpha, beta, wn, we, wd'''

        # Relabel Inputs
        pn, pe, pd, u, v, w, phi, th, psi, p, q, r = x
        delta_e, delta_a, delta_r, delta_t = delta
        w_ns, w_es, w_ds, u_wg, v_wg, w_wg = wind # w_n --> steady wind, u_w$ --> gust along $ axis

        # Compute wind data in NED
        R_bv = [[np.cos(th)*np.cos(psi), np.cos(th)*np.sin(psi), -np.sin(th)],
                [np.sin(phi)*np.sin(th)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(th)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(th)],
                [np.cos(phi)*np.sin(th)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(th)]] # R from vehicle to body frame
        V_wb = np.add(np.matmul(R_bv,[w_ns,w_es,w_ds]),[u_wg,v_wg,w_wg])
        w_n = V_wb[0]
        w_e = V_wb[1]
        w_d = V_wb[2]

        # Compute Air Data
        V_ab = [u - w_n,
                v - w_e,
                w - w_d]
        V_a = np.linalg.norm(V_ab)
        alpha = np.arctan(V_ab[2]/V_ab[0])
        beta = np.arcsin(V_ab[1]/V_a)
        self.Va = V_a
        self.alpha = alpha
        self.beta = beta

        # Compute external forces and torques on aircraft
        C_X = self.calc_Cx(alpha)
        C_Xq = self.calc_Cx_q(alpha)
        C_Xdeltae = self.calc_Cx_deltae(alpha)
        C_Z = self.calc_Cz(alpha)
        C_Zq = self.calc_Cz_q(alpha)
        C_Zdeltae = self.calc_Cz_deltae(alpha)

        prefix = 0.5 * self.rho * V_a**2 * self.S
        
        fx = -self.mass * self.g * np.sin(th) + (prefix) * (C_X + C_Xq * self.c / (2 * V_a) * q + C_Xdeltae * delta_e) + (0.5 * self.rho * self.S_prop * self.C_prop) * ((self.k_motor * delta_t)**2 - V_a**2)
        fy = self.mass * self.g * np.cos(th) * np.sin(phi) + (prefix) * (self.C_Y0 + self.C_Ybeta * beta + self.C_Yp * self.b / (2 * V_a) * p + self.C_Yr * self.b / (2 * V_a) * r + self.C_Ydeltaa * delta_a + self.C_Ydeltar * delta_r) + (0.5 * self.rho * self.S_prop * self.C_prop) * 0
        fz = self.mass * self.g * np.cos(th) * np.cos(phi) + (prefix) * (C_Z + C_Zq * self.c / (2 * V_a) * q + C_Zdeltae * delta_e) + (0.5 * self.rho * self.S_prop * self.C_prop) * 0

        l = prefix * self.b * (self.C_l0 + self.C_lbeta * beta + self.C_lp * self.b / (2 * V_a) * p + self.C_lr * self.b / (2 * V_a) * r + self.C_ldeltaa * delta_a + self.C_ldeltar * delta_r) + (-self.k_Tp * (self.k_Omega * delta_t)**2)
        m = prefix * self.c * (self.C_m0 + self.C_malpha * alpha + self.C_mq * self.c / (2 * V_a) * q + self.C_mdeltae * delta_e)
        n = prefix * self.b * (self.C_n0 + self.C_nbeta * beta + self.C_np * self.b / (2 * V_a) * p + self.C_nr * self.b / (2 * V_a) * r + self.C_ndeltaa * delta_a + self.C_ndeltar * delta_r)

        return fx,fy,fz,l,m,n,V_a,alpha,beta,w_n,w_e,w_d


    def calc_CL(self,alpha):
        sigma = (1 + np.exp(-self.M * (alpha - self.alpha_0)) + np.exp(self.M * (alpha + self.alpha_0)))/((1 + np.exp(-self.M * (alpha - self.alpha_0))) * (1 + np.exp(self.M * (alpha + self.alpha_0))))
        C_L = (1 - sigma) * (self.C_L0 + self.C_Lalpha) + (sigma * 2 * np.sign(alpha) * np.sin(alpha)**2 * np.cos(alpha))
        return C_L

    def calc_CD(self,alpha):
        AR = self.b**2/self.S
        C_D = self.C_Dp + (self.C_L0 + self.C_Lalpha * alpha)**2 / (np.pi * self.e * AR)
        return C_D

    def calc_Cx(self,alpha):
        C_D = self.calc_CD(alpha)
        C_L = self.calc_CL(alpha)
        C_x = -C_D * np.cos(alpha) + C_L * np.sin(alpha)
        return C_x

    def calc_Cx_q(self,alpha):
        Cx_q = -self.C_Dq * np.cos(alpha) + self.C_Lq * np.sin(alpha)
        return Cx_q

    def calc_Cx_deltae(self,alpha):
        Cx_deltae = -self.C_Ddeltae * np.cos(alpha) + self.C_Ldeltae * np.sin(alpha)
        return Cx_deltae

    def calc_Cz(self,alpha):
        C_D = self.calc_CD(alpha)
        C_L = self.calc_CL(alpha)
        C_Z = -C_D * np.sin(alpha) - C_L * np.cos(alpha)
        return C_Z

    def calc_Cz_q(self,alpha):
        Cz_q = -self.C_Dq * np.sin(alpha) - self.C_Lq * np.cos(alpha)
        return Cz_q

    def calc_Cz_deltae(self,alpha):
        Cz_deltae = -self.C_Ddeltae * np.sin(alpha) - self.C_Ldeltae * np.cos(alpha)
        return Cz_deltae
    
        
    def calc_dryden_gust(self,dt,description='low_low'):
        # Gust Value
        freq = np.random.randn(1)
        
        # Gust Parameters
        L_u = self.gust[description][1]
        L_v = L_u
        L_w = self.gust[description][2]
        sig_u = self.gust[description][3]
        sig_v = sig_u
        sig_w = self.gust[description][4]

        # Transfer Function numerator and denominators
        H_u_num = [sig_u * np.sqrt(2 * self.Va)]
        H_u_den = [np.sqrt(L_u), np.sqrt(L_u) * self.Va/L_u]
        
        H_v_num = [sig_v * np.sqrt(3 * self.Va),
                   sig_v * np.sqrt(3 * self.Va) * (self.Va/(np.sqrt(3) * L_v))]
        H_v_den = [np.sqrt(L_v),
                   np.sqrt(L_v) * 2 * self.Va / L_v,
                   np.sqrt(L_v) * self.Va**2/L_v**2]

        H_w_num = [sig_w * np.sqrt(3 * self.Va),
                   sig_w * np.sqrt(3 * self.Va) * (self.Va/(np.sqrt(3) * L_w))]
        H_w_den = [np.sqrt(L_w),
                   np.sqrt(L_w) * 2 * self.Va / L_w,
                   np.sqrt(L_w) * self.Va**2/L_w**2]            
    
        # Transfer Functions
        H_u = sigs.TransferFunction(H_u_num,H_u_den)
        H_v = sigs.TransferFunction(H_v_num,H_v_den)
        H_w = sigs.TransferFunction(H_w_num,H_w_den)

        # Frequency Response
        H_u = H_u.to_ss()
        H_u = sigs.lti(H_u.A,H_u.B,H_u.C,H_u.D)
        H_v = H_v.to_ss()
        H_v = sigs.lti(H_v.A,H_v.B,H_v.C,H_v.D)
        H_w = H_w.to_ss()
        H_w = sigs.lti(H_w.A,H_w.B,H_w.C,H_w.D)

        t = np.linspace(0,dt,2)
        u = [freq.tolist()[0]] * 2
        tout,u_wg,xu_wg = sigs.lsim(H_u,u,t)
        tout,v_wg,xv_wg = sigs.lsim(H_v,u,t)
        tout,w_wg,xw_wg = sigs.lsim(H_w,u,t)
        
        return u_wg[-1],v_wg[-1],w_wg[-1]
        
