from mavTrim import *

class autoPilot(Trim):
    def __init__(self):
        Trim.__init__(self)
        # Initialize Weird Variables from mavTrim that weren't initializing....
        self.variables()

        # initial Speed
        self.Va_0 = 20.0
        
        # Max Parameters
        self.max_deltas = [45*np.pi/180.0, 45*np.pi/180.0, 20*np.pi/180.0, 1.0]
        self.chi_max = 25 * np.pi / 180.0
        self.phi_max = 45 * np.pi / 180.0
        self.beta_max = 30 * np.pi / 180.0
        self.th_max = 40 * np.pi / 180.0

        # Design Parameters (I.E. These are the ones that you should tune)
        self.k_i_phi = 0.5 # roll
        self.zeta_phi = 0.9 # roll
        
        self.W_chi = 9.0 # course
        self.zeta_chi = 0.7 # course
        
        self.vr_max = 3.0 # sideslip
        self.zeta_beta = 0.7 # sideslip
        
        self.zeta_th = 0.7 # pitch

        self.W_h = 40.0 # altitude
        self.zeta_h = 0.9 # altitude

        self.W_v2 = 10.0 # airspeed pitch
        self.zeta_v2 = 1.0 # airspeed pitch

        self.zeta_v = 2.0 # airspeed throttle
        self.omega_n_v = 5.0 # airspeed throttle

        # PID Parameters
        self.tau = 5.0
        
        self.I_roll = 0.0
        self.D_roll = 0.0
        self.E_roll = 0.0

        self.I_course = 0.0
        self.E_course = 0.0
        
        self.I_side = 0.0
        self.E_side = 0.0
        
        self.D_pitch = 0.0
        self.E_pitch = 0.0

        self.I_alt = 0.0
        self.E_alt = 0.0

        self.I_air_pitch = 0.0
        self.E_air_pitch = 0.0

        self.I_air_throt = 0.0
        self.E_air_throt = 0.0

        # Initial Parameters (needed for successive loop closure)
        self.omega_n_phi = 6.0
        self.K_dc_th = 1.0
        self.omega_n_th = 10.0

        # Longitudinal State
        self.switch = 'take off'

    def sat(self,u,lim):
        # Compute Saturation
        if u > lim:
            out = lim
        elif u < -lim:
            out = -lim
        else:
            out = u
        return out

        
    def roll_attitude(self,x,phi_c,Ts):
        # Compute Transfer Function (to get coefficients)
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

        # Parameters
        self.k_p_phi = self.max_deltas[2] / self.phi_max        
        self.omega_n_phi = np.sqrt(self.k_p_phi * self.a_phi2)
        self.k_d_phi = (2 * self.zeta_phi * self.omega_n_phi - self.a_phi1) / self.a_phi2

        # PID Control Loop
        error = phi_c - phi
        if abs(error) < self.max_deltas[1]:
            self.I_roll += Ts / 2 * (error + self.E_roll)
        else:
            self.I_roll = 0
         
        self.D_roll = (2 * self.tau - Ts) / (2 * self.tau + Ts) * self.D_roll + 2 / (2 * self.tau + Ts) * (error - self.E_roll)
        # self.D_roll = -p
        self.E_roll = error

        # Get command (check for saturation)
        u = self.k_p_phi * self.E_roll + self.k_i_phi * self.I_roll + self.k_d_phi * self.D_roll
        u = self.sat(u,self.max_deltas[1])

        # Integrator Anti-Windup
        if self.k_i_phi != 0:
            u_unsat = self.k_p_phi * self.E_roll + self.k_i_phi * self.I_roll + self.k_d_phi * self.D_roll
            self.I_roll += Ts / self.k_i_phi * (u - u_unsat)

        # u is delta_a
        self.deltas[1] = u

    def course_hold(self,x,chi_c,Ts):
        # Parameters
        self.omega_n_chi = 1 / self.W_chi * self.omega_n_phi
        self.k_p_chi = 2 * self.zeta_chi * self.omega_n_chi * self.Va_0 / self.g
        self.k_i_chi = self.omega_n_chi**2 * self.Va_0 / self.g
        
        # PID Control Loop
        error = chi_c - self.chi
        if abs(error) < self.phi_max:
            self.I_course += Ts / 2 * (error + self.E_course)
        else:
            self.I_course = 0.0
        self.E_course = error

        # Get command (check for saturation)
        u = self.k_p_chi * self.E_course + self.k_i_chi * self.I_course
        u = self.sat(u,self.phi_max)
        
        # Integrator Anti-Windup
        if self.k_i_chi != 0:
            u_unsat = self.k_p_chi * self.E_course + self.k_i_chi * self.I_course
            self.I_course += Ts / self.k_i_chi * (u - u_unsat)

        # u is phi_c
        return u

    def sideslip_hold(self,x,beta_c,Ts):
        # Parameters
        # self.k_p_beta = self.max_deltas[2] / self.e_beta_max * np.sign(self.a_beta2)
        self.k_p_beta = self.max_deltas[2] / self.vr_max
        self.k_i_beta = 1 / self.a_beta2 * ((self.a_beta1 + self.a_beta2 * self.k_p_beta) / (2 * self.zeta_beta))**2
        
        # PID Control Loop
        error = beta_c - self.beta
        self.I_side += Ts / 2 * (error + self.E_side)
        self.E_side = error
        
        # Get command (check for saturation)
        u = -self.k_p_beta * self.E_side + -self.k_i_beta * self.I_side
        u = self.sat(u,self.beta_max)

        # Integrator Anti-Windup
        if self.k_i_beta != 0:
            u_unsat = -self.k_p_beta * self.E_side + -self.k_i_beta * self.I_side
            self.I_side += Ts / self.k_i_beta * (u - u_unsat)
            
        # u is delta_r
        self.deltas[2] = u
        
    def pitch_hold(self,x,th_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        
        # Parameters
        self.k_p_th = -self.max_deltas[0] / self.th_max
        self.omega_n_th = np.sqrt(self.a_th2 + self.k_p_th*self.a_th3)
        self.k_d_th = (2 * self.zeta_th * self.omega_n_th - self.a_th1) / self.a_th3
        self.K_dc_th = self.k_p_th * self.a_th3 / (self.a_th2 + self.k_p_th * self.a_th3)

        # PID Control Loop
        error = th_c - th
        self.D_pitch = (2 * self.tau - Ts) / (2 * self.tau + Ts) * self.D_pitch + 2 / (2 * self.tau + Ts) * (error - self.E_pitch)
        self.E_pitch = error

        u = self.k_p_th * self.E_pitch - self.k_d_th * q
        u = self.sat(u,self.max_deltas[0])

        # u is delta_e
        self.deltas[0] = u
        
    def altitude_hold(self,x,h_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        h = -pd

        # Parameters
        self.omega_n_h = 1 / self.W_h * self.omega_n_th
        self.k_p_h = 2 * self.zeta_h * self.omega_n_h / (self.K_dc_th * self.Va_0)
        self.k_i_h = self.omega_n_h**2 / (self.K_dc_th * self.Va_0)

        # PID Control Loop
        error = h_c - h
        if abs(error) < 25.0:
            self.I_alt += Ts / 2 * (error + self.E_alt)
        else:
            self.I_alt = 0
        self.E_alt = error

        # Get command (check for saturation)
        u = self.k_p_h * self.E_alt + self.k_i_h * self.I_alt
        u = self.sat(u,self.th_max)
        
        # Integrator Anti-Windup
        if self.k_i_h != 0:
            u_unsat = self.k_p_h * self.E_alt + self.k_i_h * self.I_alt
            self.I_alt += Ts / self.k_i_h * (u - u_unsat)

        # u is th_c
        return u

    def airspeed_pitch(self,x,v_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

        # Parameters
        self.omega_n_v2 = self.omega_n_th / self.W_v2
        self.k_i_v2 = -self.omega_n_v2**2 / (self.K_th_dc * self.g)
        self.k_p_v2 = (self.a_v1 - 2 * self.zeta_v2 * self.omega_n_v2) / (self.K_th_dc * self.g)

        # PID Control Loop
        error = v_c - self.Va
        if abs(error) < 1.0:
            self.I_air_pitch += Ts / 2 * (error + self.E_air_pitch)
        else:
            self.I_air_pitch = 0
        self.E_air_pitch = error

        # Get command (check saturation)
        u = self.k_p_v2 * self.E_air_pitch + self.k_i_v2 * self.I_air_pitch
        u = self.sat(u,self.th_max)

        if self.k_i_v2 != 0:
            u_unsat = self.k_p_v2 * self.E_air_pitch + self.k_i_v2 * self.I_air_pitch
            self.I_air_pitch += Ts / self.k_i_v2 * (u - u_unsat)

        # u is th_c
        return u

    def airspeed_throttle(self,x,v_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

        # Parameters
        self.k_p_v = (2 * self.zeta_v * self.omega_n_v - self.a_v1) / self.a_v2
        self.k_i_v = self.omega_n_v**2 / self.a_v2

        # PID Control Loop
        error = v_c - self.Va
        self.I_air_throt += Ts / 2 * (error + self.E_air_throt)
        # if abs(error) < 5.0:
        #     self.I_air_throt += Ts / 2 * (error + self.E_air_throt)
        # else:
        #     self.I_air_throt = 0
        self.E_air_throt = error

        # Get command (check saturation)
        u = self.k_p_v * self.E_air_throt + self.k_i_v * self.I_air_throt + self.delta_trim[3]
        if u > self.max_deltas[3]:
            u = self.max_deltas[3]
        elif u < 0:
            u = 0.0

        if self.k_i_v != 0:
            u_unsat = self.k_p_v * self.E_air_throt + self.k_i_v * self.I_air_throt + self.delta_trim[3]
            self.I_air_throt += Ts / self.k_i_v * (u - u_unsat)

        # u is delta_t
        self.deltas[3] = u
