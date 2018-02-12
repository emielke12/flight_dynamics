from mavTrim import *

class autoPilot(Trim):
    def __init__(self):
        Trim.__init__(self)
        
        # Max Parameters
        self.max_deltas = [45*np.pi/180.0, 45*np.pi/180.0, 45*np.pi/180.0, 1]
        self.e_th_max = 10 * np.pi / 180.0
        self.chi_max = 25 * np.pi / 180.0
        self.beta_max = 15 * np.pi / 180.0

        # Design Parameters
        self.e_phi_max = 15 * np.pi/ 180.0 # Max Roll Error
        self.k_i_phi = 0.5
        self.zeta_phi = 0.8
        self.W_chi = 20.0
        self.zeta_chi = 15.0
        self.e_beta_max = 5.0 * np.pi / 180.0 # Max Sideslip Error
        self.zeta_beta = 5.0

        # PID Parameters
        self.I_roll = 0.0
        self.D_roll = 0.0
        self.E_roll = 0.0
        self.tau_roll = 0.03
        self.I_course = 0.0
        self.E_course = 0.0
        self.I_side = 0.0
        self.E_side = 0.0
        
    def roll_attitude(self,x,phi_c,Ts):
        # Compute Transfer Function (to get coefficients)
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        self.get_transfer_functions(x)
        
        # Parameters
        self.k_p_phi = self.max_deltas[2] / self.e_phi_max * np.sign(self.a_phi2)
        self.omega_n_phi = np.sqrt(np.abs(self.a_phi2) * self.max_deltas[1] / self.e_phi_max)
        self.k_d_phi = (2 * self.zeta_phi * self.omega_n_phi - self.a_phi1) / self.a_phi2

        # PID Control Loop
        error = phi_c - phi
        self.I_roll += Ts / 2 * (error + self.E_roll)
        self.D_roll = (2 * self.tau_roll - Ts) / (2 * self.tau_roll + Ts) * self.D_roll + 2 / (2 * self.tau_roll + Ts) * (error - self.E_roll)
        self.E_roll = error

        # Get command (check for saturation)
        u_roll = self.k_p_phi * self.E_roll + self.k_i_phi * self.I_roll + self.k_d_phi * self.D_roll
        if u_roll > self.max_deltas[1]:
            u_roll = self.max_deltas[1]
        elif u_roll < -self.max_deltas[1]:
            u_roll = -self.max_deltas[1]

        if self.k_i_phi != 0:
            u_roll_unsat = self.k_p_phi * self.E_roll + self.k_i_phi * self.I_roll + self.k_d_phi * self.D_roll
            self.I_roll += Ts / self.k_i_phi * (u_roll - u_roll_unsat)

        self.deltas[1] = u_roll

    def course_hold(self,x,chi_c,Ts):
        # Parameters
        # self.omega_n_chi = 1 / self.W_chi * self.omega_n_phi
        self.omega_n_chi = 1 / self.W_chi * 11.1658
        self.k_p_chi = 2 * self.zeta_chi * self.omega_n_chi * self.Vg / self.g
        self.k_i_chi = self.omega_n_chi**2 * self.Vg / self.g

        # PID Control Loop
        error = chi_c - self.chi
        self.I_course += Ts / 2 * (error + self.E_course)
        self.E_course = error

        # Get command (check for saturation)
        u_course = self.k_p_chi * self.E_course + self.k_i_chi * self.I_course
        if u_course > self.chi_max:
            u_course = self.chi_max
        elif u_course < -self.chi_max:
            u_course = -self.chi_max

        if self.k_i_chi != 0:
            u_course_unsat = self.k_p_chi * self.E_course + self.k_i_chi * self.I_course
            self.I_course += Ts / self.k_i_chi * (u_course - u_course_unsat)

        # u_course is phi_c
        return u_course

    def sideslip_hold(self,x,beta_c,Ts):
        # Parameters
        self.k_p_beta = self.max_deltas[2] / self.e_beta_max * np.sign(self.a_beta2)
        self.k_i_beta = 1 / self.a_beta2 * ((self.a_beta1 + self.a_beta2 * self.k_p_beta) / (2 * self.zeta_beta))**2
        
        # PID Control Loop
        error = beta_c - beta
        self.I_side += Ts / 2 * (error + self.E_side)
        self.E_side = error
        
        # Get command (check for saturation)
        u_side = self.k_p_beta * self.E_side + self.k_i_beta * self.I_side
        if u_side > self.beta_max:
            u_side = self.beta_max
        elif u_side < -self.beta_max:
            u_side = -self.beta_max

        if self.k_i_beta != 0:
            u_side_unsat = self.k_p_beta * self.E_side + self.k_i_beta * self.I_side
            self.I_side += Ts / self.k_i_beta * (u_side - u_side_unsat)
            
        # u_side is delta_r
        self.deltas[2] = u_side
        
        
        
