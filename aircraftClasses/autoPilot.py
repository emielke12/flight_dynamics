from mavTrim import *

class autoPilot(Trim):
    def __init__(self):
        Trim.__init__(self)

        self.max_deltas = [45*np.pi/180.0, 45*np.pi/180.0, 45*np.pi/180.0, 1]
        self.phi_max = 15 * np.pi / 180.0
        self.e_th_max = 10 * np.pi / 180.0

        # Design Parameters
        self.e_phi_max = 20 * np.pi/ 180.0
        self.zeta_phi = 1.0
        
    def roll_attitude(self,x,phi_c):
        # Compute Transfer Function (to get coefficients)
        self.get_transfer_function(x)

        # Parameters
        self.k_p_phi = self.max_deltas[2] / self.e_phi_max * np.sign(self.a_phi2)
        self.omega_n_phi = np.sqrt(np.abs(self.a_phi2) * self.max_deltas[1] / self.e_phi_max)
        self.k_d_phi = (2 * self.zeta_phi * self.omega_n_phi - self.a_phi1) / self.a_phi2
        self.k_i_phi = 0.01

        self.deltas[1] = 
        
        
