from kalmanFilter.py import *

class guidanceModel(kalmanFilter):
    def __init__(self,x0,trim):
        kalmanFilter.__init__(self,x0,trim)

        # Parameters to tune
        self.bchidot = 1.0
        self.bchi = 1.0
        self.bhdot = 1.0
        self.bh = 1.0
        self.bva = 1.0
        self.bphi = 1.0
        
    def model1(self,x,wind,commands):
        # Rename Inputs
        pn,pe,chi,chidot,h,hdot,va = x
        wn,we,wd = wind
        chidot_c, chi_c, hdot_c, h_c, va_c = commands

        # Eq 9.19
        pndot = self.Va_hat * np.cos(psi) + wn
        pedot = self.Va_hat * np.sin(psi) + we
        chiddot = self.bchidot * (chidot_c - chidot) + self.bchi * (chi_c - chi)
        hddot = self.bhdot * (hdot_c - hdot) + self.bh * (h_c - h)
        vadot = self.bva * (va_c - va)

        return [pndot,pedot,chidot,chiddot,hdot,hddot,vadot]

    def model2(self,x,wind,commands):
        # Rename Inputs
        pn, pe, psi, h, hdot, va, phi = x
        wn, we, wd = wind
        chidot_c, chi_c, hdot_c, h_c, va_c, phi_c = commands

        # Eq 9.20
        pndot = self.Va_hat * np.cos(psi) + wn
        pedot = self.Va_hat * np.sin(psi) + we
        psidot = self.g / self.Va_hat * np.tan(phi)
        hddot = self.bhdot * (hdot_c - hdot) + self.bh * (h_c - h)
        vadot = self.bva * (va_c - va)
        phidot = self.bphi * (phi_c - phi)

        return [pndot,pedot,psidot,hdot,hddot,vadot,phidot]
