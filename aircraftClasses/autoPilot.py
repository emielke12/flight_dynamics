from mavTrim import *

class autoPilot(Trim):
    def __init__(self,x0,trim):
        Trim.__init__(self)
        # Initialize Weird Variables from mavTrim that weren't initializing....
        self.variables()

        # Set Trim values
        self.trim = trim
        
        # Get Trim values
        x0 = self.call_opt(x0,[0.0,0.0,0.0,0.0,0.0,0.0])
        
        # Transfer functions
        self.get_transfer_functions(x0)
        
        # initial things
        self.Va_0 = np.linalg.norm([x0[3:6]])
        self.h_takeoff = 15.0
        self.h_hold = 60.0
        self.h_band = 20.0
        self.h_zone = 0.5
        self.th_takeoff = 25.0 * np.pi / 180.0
        self.switch = 'takeoff'
        
        # Max Parameters
        self.max_deltas = [45*np.pi/180.0, 45*np.pi/180.0, 20*np.pi/180.0, 1.0]
        self.chi_max = 25 * np.pi / 180.0
        self.phi_max = 45 * np.pi / 180.0
        self.beta_max = 30 * np.pi / 180.0
        self.th_max = 50 * np.pi / 180.0

        # Design Parameters (I.E. These are the ones that you should tune)
        # Roll
        self.k_i_phi = 0.5 
        self.zeta_phi = 1.8 
        self.k_p_phi = self.max_deltas[2] / self.phi_max 
        self.omega_n_phi = np.sqrt(self.k_p_phi * self.a_phi2)
        self.k_d_phi = (2 * self.zeta_phi * self.omega_n_phi - self.a_phi1) / self.a_phi2

        # Course
        self.W_chi = 17.0
        self.zeta_chi = 1.5
        self.omega_n_chi = 1 / self.W_chi * self.omega_n_phi
        self.k_p_chi = 2 * self.zeta_chi * self.omega_n_chi * self.Va_0 / self.g
        self.k_i_chi = self.omega_n_chi**2 * self.Va_0 / self.g

        # Sideslip
        self.vr_max = 3.0 
        self.zeta_beta = 0.7
        # self.k_p_beta = self.max_deltas[2] / self.e_beta_max * np.sign(self.a_beta2)
        self.k_p_beta = self.max_deltas[2] / self.vr_max
        self.k_i_beta = 1 / self.a_beta2 * ((self.a_beta1 + self.a_beta2 * self.k_p_beta) / (2 * self.zeta_beta))**2

        # Pitch
        self.zeta_th = 0.7
        self.k_p_th = -self.max_deltas[0] / self.th_max * 1.0
        self.omega_n_th = np.sqrt(self.a_th2 + self.k_p_th*self.a_th3)
        self.k_d_th = (2 * self.zeta_th * self.omega_n_th - self.a_th1) / self.a_th3
        self.K_dc_th = self.k_p_th * self.a_th3 / (self.a_th2 + self.k_p_th * self.a_th3)

        # Altitude
        self.W_h = 25.0 
        self.zeta_h = 1.5
        self.omega_n_h = 1 / self.W_h * self.omega_n_th
        self.k_p_h = 2 * self.zeta_h * self.omega_n_h / (self.K_dc_th * self.Va_0)
        self.k_i_h = self.omega_n_h**2 / (self.K_dc_th * self.Va_0)

        # Airspeed Pitch
        self.W_v2 = 8.0 
        self.zeta_v2 = 1.0
        self.omega_n_v2 = self.omega_n_th / self.W_v2
        self.k_i_v2 = -self.omega_n_v2**2 / (self.K_dc_th * self.g)
        self.k_p_v2 = (self.a_v1 - 2 * self.zeta_v2 * self.omega_n_v2) / (self.K_dc_th * self.g)

        # Airspeed Throttle
        self.zeta_v = 0.7 
        self.omega_n_v = 5.0
        self.k_p_v = (2 * self.zeta_v * self.omega_n_v - self.a_v1) / self.a_v2
        self.k_i_v = self.omega_n_v**2 / self.a_v2

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

        self.I_alt = 30.5
        self.E_alt = 0.0

        self.I_air_pitch = 0.0
        self.E_air_pitch = 0.0

        self.I_air_throt = 0.0
        self.E_air_throt = 0.0

    def sat(self,u,ulim,llim):
        # Compute Saturation
        if u > ulim:
            out = ulim
        elif u < llim:
            out = llim
        else:
            out = u
        return out

        
    def roll_attitude(self,x,phi_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

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
        u = self.sat(u,self.max_deltas[1],-self.max_deltas[1])

        # Integrator Anti-Windup
        if self.k_i_phi != 0:
            u_unsat = self.k_p_phi * self.E_roll + self.k_i_phi * self.I_roll + self.k_d_phi * self.D_roll
            self.I_roll += Ts / self.k_i_phi * (u - u_unsat)

        # u is delta_a
        self.deltas[1] = u

    def course_hold(self,x,chi_c,Ts):        
        # PID Control Loop
        error = chi_c - self.chi
        if abs(error) < self.chi_max:
            self.I_course += Ts / 2 * (error + self.E_course)
        else:
            self.I_course = 0.0
        self.E_course = error

        # Get command (check for saturation)
        u = self.k_p_chi * self.E_course + self.k_i_chi * self.I_course
        u = self.sat(u,self.phi_max,-self.phi_max)
        
        # Integrator Anti-Windup
        if self.k_i_chi != 0:
            u_unsat = self.k_p_chi * self.E_course + self.k_i_chi * self.I_course
            self.I_course += Ts / self.k_i_chi * (u - u_unsat)

        # u is phi_c
        self.roll_attitude(x,u,Ts)

    def sideslip_hold(self,x,Ts):
        # PID Control Loop
        error = -self.beta
        self.I_side += Ts / 2 * (error + self.E_side)
        self.E_side = error
        
        # Get command (check for saturation)
        u = -self.k_p_beta * self.E_side + -self.k_i_beta * self.I_side
        u = self.sat(u,self.beta_max,-self.beta_max)

        # Integrator Anti-Windup
        if self.k_i_beta != 0:
            u_unsat = -self.k_p_beta * self.E_side + -self.k_i_beta * self.I_side
            self.I_side += Ts / self.k_i_beta * (u - u_unsat)
            
        # u is delta_r
        self.deltas[2] = u
        
    def pitch_hold(self,x,th_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        
        # PID Control Loop
        error = th_c - th
        self.D_pitch = (2 * self.tau - Ts) / (2 * self.tau + Ts) * self.D_pitch + 2 / (2 * self.tau + Ts) * (error - self.E_pitch)
        self.E_pitch = error

        u = self.k_p_th * self.E_pitch - self.k_d_th * q
        u = self.sat(u,self.max_deltas[0],-self.max_deltas[0])

        # u is delta_e
        self.deltas[0] = u
        
    def altitude_hold(self,x,h_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        h = -pd

        # PID Control Loop
        error = h_c - h
        if abs(error) < 20.0:
            self.I_alt += Ts / 2 * (error + self.E_alt)
        else:
            self.I_alt = 0
        self.E_alt = error

        # Get command (check for saturation)
        u = self.k_p_h * self.E_alt + self.k_i_h * self.I_alt
        u = self.sat(u,self.th_max,-self.th_max)
        
        # Integrator Anti-Windup
        if self.k_i_h != 0:
            u_unsat = self.k_p_h * self.E_alt + self.k_i_h * self.I_alt
            self.I_alt += Ts / self.k_i_h * (u - u_unsat)

        # u is th_c
        self.pitch_hold(x,u,Ts)
        # return u

    def airspeed_pitch(self,x,v_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

        # PID Control Loop
        error = v_c - self.Va
        if abs(error) < 3.0:
            self.I_air_pitch += Ts / 2 * (error + self.E_air_pitch)
        else:
            self.I_air_pitch = 0
        self.E_air_pitch = error

        # Get command (check saturation)
        u = self.k_p_v2 * self.E_air_pitch + self.k_i_v2 * self.I_air_pitch
        u = self.sat(u,self.th_max,-self.th_max)

        # Integrator Anti-Windup
        if self.k_i_v2 != 0:
            u_unsat = self.k_p_v2 * self.E_air_pitch + self.k_i_v2 * self.I_air_pitch
            self.I_air_pitch += Ts / self.k_i_v2 * (u - u_unsat)

        # u is th_c
        self.pitch_hold(x,u,Ts)
        # return u

    def airspeed_throttle(self,x,v_c,Ts):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

        # PID Control Loop
        error = v_c - self.Va
        self.I_air_throt += Ts / 2 * (error + self.E_air_throt)
        if abs(error) < 1.5:
            self.I_air_throt += Ts / 2 * (error + self.E_air_throt)
        else:
            self.I_air_throt = 0
        self.E_air_throt = error

        # Get command (check saturation)
        u = self.k_p_v * self.E_air_throt + self.k_i_v * self.I_air_throt + self.delta_trim[3]
        u = self.sat(u,self.max_deltas[3],0)


        # Integrator Anti-Windup
        if self.k_i_v != 0:
            u_unsat = self.k_p_v * self.E_air_throt + self.k_i_v * self.I_air_throt + self.delta_trim[3]
            self.I_air_throt += Ts / self.k_i_v * (u - u_unsat)

        # u is delta_t
        self.deltas[3] = u

    def mode(self,x):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        h = -pd
        
        if self.switch == 'takeoff':
            print 'Takeoff'
            if h > self.h_takeoff:
                self.switch = 'air'
            
            self.deltas[3] = 1.0
            th_c = deepcopy(self.th_takeoff)
            self.pitch_hold(x,th_c,self.dt)

        if self.switch == 'air':
            # Set Throttle
            if h < self.h_hold - self.h_band:
                self.deltas[3] = 1.0 # Climb
                # print 'Climbing'
            elif h > self.h_hold + self.h_band:
                self.deltas[3] = 0.0 # Descend
                # print 'Descending'

            if h > self.h_hold - self.h_band and h < self.h_hold + self.h_band:
                # print 'Hold Zone: Altitude = ',self.h_hold#,'+/-',self.h_band,'meters'
                # Altitude Hold with Throttle
                self.airspeed_throttle(x,self.Va_0,self.dt)
                self.altitude_hold(x,self.h_hold,self.dt)
                # self.pitch_hold(x,th_c,self.dt)
            else:
                # Climb and Descend
                self.airspeed_pitch(x,self.Va_0,self.dt)
                # self.pitch_hold(x,th_c,self.dt)

            
