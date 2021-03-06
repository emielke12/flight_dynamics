from sensors import *
from copy import deepcopy
import time, sys
import pyqtgraph as pg

class kalmanFilter(Sensors):
    def __init__(self,x0,trim):
        Sensors.__init__(self,x0,trim)
        self.sens = [0.0] * 8
        self.gps = [0.0] * 5
        self.dt = 0.01
        self.N = 10
        self.Tout = self.dt
        np.set_printoptions(precision = 2, suppress = True)
        self.x_hat = x0
        
        # Sensed States
        self.p_hat = 0.0
        self.q_hat = 0.0
        self.r_hat = 0.0
        self.h_hat = -x0[2]
        self.stat_hat = -x0[2] * self.rho * self.g
        self.Va_hat = np.sqrt(x0[3]**2 + x0[4]**2 + x0[5]**2)
        self.diff_hat = self.Va_hat**2 * self.rho / 2
        self.phi_hat = 0.0
        self.th_hat = 0.0
        self.pn_hat = 0.0
        self.pe_hat = 0.0
        self.chi_hat = 0.0
        self.psi_hat = 0.0
        self.Vg_hat = np.linalg.norm(x0[3:6])

        # Mid States
        self.stat_hat = -x0[2]
        self.diff_hat = 0.0
        self.ax_hat = 0.0
        self.ay_hat = 0.0
        self.az_hat = 0.0
        
        # Low Pass Filter alpha values
        self.alpha_p = 0.5
        self.alpha_q = 0.5
        self.alpha_r = 0.5
        self.alpha_h = 0.5
        self.alpha_va = 0.5
        self.alpha_a = 0.8
        self.alpha_pn = 0.5
        self.alpha_pe = 0.5
        self.alpha_chi = 0.5
        self.alpha_vg = 0.5

        # Previous Sensor Values
        self.y_attprev = [0.0, 0.0, 0.0]
        self.y_gpsprev = [0.0, 0.0, 0.0, 0.0]
        
        # Covariances and Kalman Filter Stuff
        self.x_atthat = [0.0, 0.0] # phi, theta
        self.Pa = np.diag([(15.0*np.pi/180.0)**2, (15.0*np.pi/180.0)**2])
        self.Qa = np.diag([1e-9,1e-9])
        self.Ra = np.diag([1e6,1e6,1e6])
        
        self.x_gpshat = [0.0, 0.0, self.Vg_hat, 0.0, 0.0] # pn,pe,vg,chi,psi
        self.Pg = np.diag([100, 100, 1.0, (10.0 * np.pi/180.0)**2, (10.0 * np.pi/180.0)**2])
        self.Qg = np.diag([0.01, 0.01, 10.0, 10.0, 0.00001])
        self.Rg = np.diag([self.sigma_n**2,self.sigma_e**2,(self.sigma_V/10)**2,self.sigma_chi**2])

    def ekf(self,x,f,wind,count,plot=True):
        # Rename inputs
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        wn,we,wd = wind
        fx,fy,fz = f

        # Get Sensor Readings
        self.get_lpf(x,fx,fy,fz,wn,we,wd,count,False)
        
        # roll and pitch stuff
        u_att = [self.p_hat,self.q_hat,self.r_hat,self.Va_hat]
        y_att = self.sens[3:6]
        x_atthat,P_att = self.predict(self.x_atthat,self.f_att,u_att,self.df_att,self.Pa,self.Qa)
        for i in xrange(len(y_att)):
            if y_att[i] != self.y_attprev[i]:
                x_atthat,P_att = self.update(x_atthat,self.dh_att,u_att,P_att,self.Ra,y_att[i],self.h_att,i)
                self.y_attprev[i] = y_att[i]
        
        self.Pa = P_att
        self.x_atthat = x_atthat
        self.phi_hat,self.th_hat = x_atthat

        # Gps stuff
        u_gps = [self.Va_hat, self.q_hat, self.r_hat, self.phi_hat, self.th_hat, wn, we]
        y_gps = self.gps[:2] + self.gps[3:]
        x_gpshat, P_gps = self.predict(self.x_gpshat,self.f_gps,u_gps,self.df_gps,self.Pg,self.Qg)
        for i in xrange(len(y_gps)):
            if y_gps[i] != self.y_gpsprev[i]:
                x_gpshat,P_gps = self.update(x_gpshat,self.dh_gps,u_gps,P_gps,self.Rg,y_gps[i],self.h_gps,i)
                self.y_gpsprev[i] = y_gps[i]

        self.Pg = P_gps
        self.x_gpshat = x_gpshat
        self.pn_hat = x_gpshat[0]
        self.pe_hat = x_gpshat[1]
        self.psi_hat = x_gpshat[4]
        self.Vg_hat = x_gpshat[2]
        self.chi_hat = x_gpshat[3]

        # Plot
        if plot == True:
            self.ekfplot2(np.concatenate((x_atthat,x_gpshat)),[phi,th,pn,pe,self.Vg,self.chi,psi],count)

        # Set to hat variables
        self.x_hat = [self.pn_hat, self.pe_hat, -self.h_hat, u, v, w, self.phi_hat, self.th_hat, self.psi_hat,self.p_hat, self.q_hat, self.r_hat] 

    def predict(self,xhat,f,u,df,P,Q):
        # EKF Prediction Step
        for i in xrange(self.N):
            xhat += np.multiply(self.Tout / self.N, f(xhat,u))
            A = df(xhat,u)
            P += self.Tout / self.N * (np.add(np.add(np.matmul(A, P), np.matmul(P, np.transpose(A))), Q))
        return xhat,P

    def update(self,xhat,dh,u,P,R,y,h,row):
        # EKF Measurment Update
        C = dh(xhat,u,row)
        L = np.matmul(P,np.transpose(C)) / (R[row][row] + np.matmul(np.matmul(C,P),np.transpose(C)))
        P = (np.eye(np.shape(L * C)[0]) - L * C) * P
        xhat += L * (y - h(xhat,u)[row])
        return xhat,P

    # Attitude f, df/dx, h, dh/dx models
    def f_att(self,x,u):
        phi, th = x
        p, q, r, va = u
        f = [p + q * np.sin(phi) * np.tan(th) + r * np.cos(phi) * np.tan(th),
             q * np.cos(phi) - r * np.sin(phi)]
        return f

    def h_att(self,x,u):
        phi, th = x
        p, q, r, va = u
        f = [q * va * np.sin(th) + self.g * np.sin(th),
             r * va * np.cos(th) - p * va * np.sin(th) - self.g * np.cos(th) * np.sin(phi),
             -q * va * np.cos(th) - self.g * np.cos(th) * np.cos(phi)]
        return f

    def df_att(self,x,u):
        phi, th = x
        p, q, r, va = u
        f = [[q * np.cos(phi) * np.tan(th) - r * np.sin(phi) * np.tan(th), (q * np.sin(phi) + r * np.cos(phi)) / (np.cos(th)**2)],
             [-q * np.sin(phi) - r * np.cos(phi), 0]]
        return f

    def dh_att(self,x,u,row):
        phi, th = x
        p, q, r, va = u
        f = [[0, q * va * np.cos(th) + self.g * np.cos(th)],
             [-self.g * np.cos(phi) * np.cos(th), -r * va * np.sin(th) - p * va * np.cos(th) + self.g * np.sin(phi) * np.sin(th)],
             [self.g * np.sin(phi) * np.cos(th), (q * va + self.g * np.cos(phi)) * np.sin(th)]]
        return f[row]

    # GPS f, df/dx, h, dh/dx models
    def f_gps(self,x,u):
        pn, pe, vg, chi, psi= x
        va, q, r, phi, th, wn, we = u
        psidot = q * np.sin(phi) / np.cos(th) + r * np.cos(phi) / np.cos(th)
        f = [vg * np.cos(chi),
             vg * np.sin(chi),
             psidot * va * np.sin(chi - psi),
             self.g / vg * np.tan(phi) * np.cos(chi - psi),
             psidot]
        return f

    def h_gps(self,x,u):
        pn, pe, vg, chi, psi = x
        va, q, r, phi, th, _, _ = u
        f = [pn,
             pe,
             vg,
             chi]
        return f
    
    def df_gps(self,x,u):
        pn, pe, vg, chi, psi = x
        va, q, r, phi, th, wn, we = u
        psidot = q * (np.sin(phi) / np.cos(th)) + r * (np.cos(phi) / np.cos(th))
        # vgdot = va / vg * psidot * (-wn * np.sin(psi) + we * np.cos(psi))
        vgdot = psidot * va * np.sin(chi - psi)
        f = [[0, 0, np.cos(chi), -vg * np.sin(chi), 0],
             [0, 0, np.sin(chi), vg * np.cos(chi), 0],
             [0, 0, -vgdot / vg, psidot * va * np.cos(chi - psi), -psidot * va * np.cos(chi - psi)],
             [0, 0, -self.g / (vg**2) * np.tan(phi) * np.cos(chi - psi), -self.g / vg * np.tan(phi) * np.sin(chi - psi), self.g / vg * np.tan(phi) * np.sin(chi - psi)],
             [0, 0, 0, 0, 0]]
        return f

    def dh_gps(self,x,u,row):
        f = [[1, 0, 0, 0, 0],
             [0, 1, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [0, 0, 0, 1, 0]]
        return f[row]
    
    # Low-pass Filter 
    def get_lpf(self,x,fx,fy,fz,wn,we,wd,count,plot = True):
        # rename inputs
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        
        # Get Sensor Readings
        if self.t_sim % self.Ts < self.dt:
            self.sens = self.sensor_readings(x,fx,fy,fz)
        if self.t_sim % self.Ts_gps < self.dt:
            self.gps = self.gps_readings(x,wn,we,wd)
        y_gyro_x,y_gyro_y,y_gyro_z,y_accel_x,y_accel_y,y_accel_z,y_abs_pres,y_diff_pres = self.sens
        y_gps_n,y_gps_e,y_gps_d,y_gps_Vg,y_gps_chi = self.gps

        # Pass certain ones through low pass filtering
        # self.p_hat = self.lpf(self.p_hat,y_gyro_x,self.alpha_p)
        # self.q_hat = self.lpf(self.q_hat,y_gyro_y,self.alpha_q)
        # self.r_hat = self.lpf(self.r_hat,y_gyro_z,self.alpha_r)
        self.p_hat = y_gyro_x
        self.q_hat = y_gyro_y
        self.r_hat = y_gyro_z
        self.stat_hat = self.lpf(self.stat_hat,y_abs_pres,self.alpha_h)
        self.h_hat = self.stat_hat / (self.rho * self.g)
        self.diff_hat = self.lpf(self.diff_hat,y_diff_pres,self.alpha_va)
        self.Va_hat = np.sqrt(2 / self.rho * self.diff_hat)
        self.ax_hat = self.lpf(self.ax_hat, y_accel_x, self.alpha_a)
        self.ay_hat = self.lpf(self.ay_hat, y_accel_y, self.alpha_a)
        self.az_hat = self.lpf(self.az_hat, y_accel_z, self.alpha_a)
        # self.ax_hat = self.g * np.sin(th)
        # self.ay_hat = -self.g * np.cos(th) * np.sin(phi)
        # self.az_hat = -self.g * np.cos(th) * np.cos(phi)
        # self.phi_hat = np.arctan(self.ay_hat / self.az_hat)
        # self.th_hat = np.arcsin(self.ax_hat / self.g)
        self.pn_hat = self.lpf(self.pn_hat, y_gps_n, self.alpha_pn)
        self.pe_hat = self.lpf(self.pe_hat, y_gps_e, self.alpha_pe)
#         self.chi_hat = self.lpf(self.chi_hat, y_gps_chi, self.alpha_chi)
#         self.Vg_hat = self.lpf(self.Vg_hat, y_gps_Vg, self.alpha_vg)

        # Plot
        s = [self.pn_hat,self.pe_hat,self.h_hat,self.Vg_hat,self.chi_hat,
             self.p_hat,self.q_hat,self.r_hat,self.Va_hat,self.phi_hat,self.th_hat]
        true = [pn,pe,-pd,self.Vg,self.chi,p,q,r,self.Va,phi,th]
        if plot == True:
            self.update_sensor_plot(s,true,count)
        
    def lpf(self,y,u,alpha):
        # Unity DC-gain Low-Pass Filter
        return alpha * y + (1 - alpha) * u

    def update_sensor_plot(self,s,true,count):
        if self.sensor_first == True:
            self.sensor_fig = plt.figure(figsize=(20,10))
            self.sax1 = self.sensor_fig.add_subplot(431)
            self.sax2 = self.sensor_fig.add_subplot(432)
            self.sax3 = self.sensor_fig.add_subplot(433)
            self.sax4 = self.sensor_fig.add_subplot(434)
            self.sax5 = self.sensor_fig.add_subplot(435)
            self.sax6 = self.sensor_fig.add_subplot(436)
            self.sax7 = self.sensor_fig.add_subplot(437)
            self.sax8 = self.sensor_fig.add_subplot(438)
            self.sax9 = self.sensor_fig.add_subplot(439)
            self.sax10 = self.sensor_fig.add_subplot(4,3,10)
            self.sax11 = self.sensor_fig.add_subplot(4,3,11)

            self.sensor_fig.show()
            self.sensor_first = False
            self.axis_font = {'fontsize':'18'}
            
            self.sax1.set_ylabel(r'$\hat{p}_{n} (m)$',**self.axis_font)
            self.sax2.set_ylabel(r'$\hat{p}_{e} (m)$',**self.axis_font)
            self.sax3.set_ylabel(r'$\hat{h} (m)$',**self.axis_font)
            self.sax4.set_ylabel(r'$\hat{V}_g (m/s)$',**self.axis_font)
            self.sax5.set_ylabel(r'$\hat{\chi} (rad)$',**self.axis_font)
            self.sax6.set_ylabel(r'$\hat{p} (rad/s)$',**self.axis_font)
            self.sax7.set_ylabel(r'$\hat{q} (rad/s)$',**self.axis_font)
            self.sax8.set_ylabel(r'$\hat{r} (rad/s)$',**self.axis_font)
            self.sax9.set_ylabel(r'$\hat{V}_a (m/s)$',**self.axis_font)
            self.sax10.set_ylabel(r'$\hat{\phi} (rad)$',**self.axis_font)
            self.sax11.set_ylabel(r'$\hat{\theta} (rad)$',**self.axis_font)

            self.s1, = self.sax1.plot(s[0])
            self.s2, = self.sax2.plot(s[1])
            self.s3, = self.sax3.plot(s[2])
            self.s4, = self.sax4.plot(s[3])
            self.s5, = self.sax5.plot(s[4])
            self.s6, = self.sax6.plot(s[5])
            self.s7, = self.sax7.plot(s[6])
            self.s8, = self.sax8.plot(s[7])
            self.s9, = self.sax9.plot(s[8])
            self.s10, = self.sax10.plot(s[9])
            self.s11, = self.sax11.plot(s[10])
            self.t1, = self.sax1.plot(true[0])
            self.t2, = self.sax2.plot(true[1])
            self.t3, = self.sax3.plot(true[2])
            self.t4, = self.sax4.plot(true[3])
            self.t5, = self.sax5.plot(true[4])
            self.t6, = self.sax6.plot(true[5])
            self.t7, = self.sax7.plot(true[6])
            self.t8, = self.sax8.plot(true[7])
            self.t9, = self.sax9.plot(true[8])
            self.t10, = self.sax10.plot(true[9])
            self.t11, = self.sax11.plot(true[10])

        else:
            t = np.linspace(0,self.t_sim,count)

            self.s1.set_xdata(t)
            self.s2.set_xdata(t)
            self.s3.set_xdata(t)
            self.s4.set_xdata(t)
            self.s5.set_xdata(t)
            self.s6.set_xdata(t)
            self.s7.set_xdata(t)
            self.s8.set_xdata(t)
            self.s9.set_xdata(t)
            self.s10.set_xdata(t)
            self.s11.set_xdata(t)
            self.t1.set_xdata(t)
            self.t2.set_xdata(t)
            self.t3.set_xdata(t)
            self.t4.set_xdata(t)
            self.t5.set_xdata(t)
            self.t6.set_xdata(t)
            self.t7.set_xdata(t)
            self.t8.set_xdata(t)
            self.t9.set_xdata(t)
            self.t10.set_xdata(t)
            self.t11.set_xdata(t)

            self.s1.set_ydata(np.append(self.s1.get_ydata(),s[0]))
            self.s2.set_ydata(np.append(self.s2.get_ydata(),s[1]))
            self.s3.set_ydata(np.append(self.s3.get_ydata(),s[2]))
            self.s4.set_ydata(np.append(self.s4.get_ydata(),s[3]))
            self.s5.set_ydata(np.append(self.s5.get_ydata(),s[4]))
            self.s6.set_ydata(np.append(self.s6.get_ydata(),s[5]))
            self.s7.set_ydata(np.append(self.s7.get_ydata(),s[6]))
            self.s8.set_ydata(np.append(self.s8.get_ydata(),s[7]))
            self.s9.set_ydata(np.append(self.s9.get_ydata(),s[8]))
            self.s10.set_ydata(np.append(self.s10.get_ydata(),s[9]))
            self.s11.set_ydata(np.append(self.s11.get_ydata(),s[10]))
            self.t1.set_ydata(np.append(self.t1.get_ydata(),true[0]))
            self.t2.set_ydata(np.append(self.t2.get_ydata(),true[1]))
            self.t3.set_ydata(np.append(self.t3.get_ydata(),true[2]))
            self.t4.set_ydata(np.append(self.t4.get_ydata(),true[3]))
            self.t5.set_ydata(np.append(self.t5.get_ydata(),true[4]))
            self.t6.set_ydata(np.append(self.t6.get_ydata(),true[5]))
            self.t7.set_ydata(np.append(self.t7.get_ydata(),true[6]))
            self.t8.set_ydata(np.append(self.t8.get_ydata(),true[7]))
            self.t9.set_ydata(np.append(self.t9.get_ydata(),true[8]))
            self.t10.set_ydata(np.append(self.t10.get_ydata(),true[9]))
            self.t11.set_ydata(np.append(self.t11.get_ydata(),true[10]))

            self.sax1.relim()
            self.sax2.relim()
            self.sax3.relim()
            self.sax4.relim()
            self.sax5.relim()
            self.sax6.relim()
            self.sax7.relim()
            self.sax8.relim()
            self.sax9.relim()
            self.sax10.relim()
            self.sax11.relim()
            self.sax1.autoscale_view()
            self.sax2.autoscale_view()
            self.sax3.autoscale_view()
            self.sax4.autoscale_view()
            self.sax5.autoscale_view()
            self.sax6.autoscale_view()
            self.sax7.autoscale_view()
            self.sax8.autoscale_view()
            self.sax9.autoscale_view()
            self.sax10.autoscale_view()
            self.sax11.autoscale_view()

            
            self.sensor_fig.tight_layout()
            self.sensor_fig.canvas.draw()

    def ekfplot(self,s,true,count):
        if self.sensor_first == True:
            self.sensor_fig = plt.figure(figsize=(20,10))
            self.sax1 = self.sensor_fig.add_subplot(331)
            self.sax2 = self.sensor_fig.add_subplot(332)
            self.sax3 = self.sensor_fig.add_subplot(333)
            self.sax4 = self.sensor_fig.add_subplot(334)
            self.sax5 = self.sensor_fig.add_subplot(335)
            self.sax6 = self.sensor_fig.add_subplot(336)
            self.sax7 = self.sensor_fig.add_subplot(337)
            self.sax8 = self.sensor_fig.add_subplot(338)
            self.sax9 = self.sensor_fig.add_subplot(339)

            self.sensor_fig.show()
            self.sensor_first = False
            self.axis_font = {'fontsize':'18'}
            
            self.sax1.set_ylabel(r'$\phi (rad)$',**self.axis_font)
            self.sax2.set_ylabel(r'$\theta (rad)$',**self.axis_font)
            self.sax3.set_ylabel(r'$p_n (m)$',**self.axis_font)
            self.sax4.set_ylabel(r'$p_e (m)$',**self.axis_font)
            self.sax5.set_ylabel(r'$V_g (m/s)$',**self.axis_font)
            self.sax6.set_ylabel(r'$\chi (rad)$',**self.axis_font)
            self.sax9.set_ylabel(r'$\psi (rad)$',**self.axis_font)

            self.s1, = self.sax1.plot(s[0])
            self.s2, = self.sax2.plot(s[1])
            self.s3, = self.sax3.plot(s[2])
            self.s4, = self.sax4.plot(s[3])
            self.s5, = self.sax5.plot(s[4])
            self.s6, = self.sax6.plot(s[5])
            self.s7, = self.sax7.plot(s[6])
            self.s8, = self.sax8.plot(s[7])
            self.s9, = self.sax9.plot(s[8])
            self.t1, = self.sax1.plot(true[0])
            self.t2, = self.sax2.plot(true[1])
            self.t3, = self.sax3.plot(true[2])
            self.t4, = self.sax4.plot(true[3])
            self.t5, = self.sax5.plot(true[4])
            self.t6, = self.sax6.plot(true[5])
            self.t7, = self.sax7.plot(true[6])
            self.t8, = self.sax8.plot(true[7])
            self.t9, = self.sax9.plot(true[8])

        else:
            t = np.linspace(0,self.t_sim,count)

            self.s1.set_xdata(t)
            self.s2.set_xdata(t)
            self.s3.set_xdata(t)
            self.s4.set_xdata(t)
            self.s5.set_xdata(t)
            self.s6.set_xdata(t)
            self.s7.set_xdata(t)
            self.s8.set_xdata(t)
            self.s9.set_xdata(t)
            self.t1.set_xdata(t)
            self.t2.set_xdata(t)
            self.t3.set_xdata(t)
            self.t4.set_xdata(t)
            self.t5.set_xdata(t)
            self.t6.set_xdata(t)
            self.t7.set_xdata(t)
            self.t8.set_xdata(t)
            self.t9.set_xdata(t)

            self.s1.set_ydata(np.append(self.s1.get_ydata(),s[0]))
            self.s2.set_ydata(np.append(self.s2.get_ydata(),s[1]))
            self.s3.set_ydata(np.append(self.s3.get_ydata(),s[2]))
            self.s4.set_ydata(np.append(self.s4.get_ydata(),s[3]))
            self.s5.set_ydata(np.append(self.s5.get_ydata(),s[4]))
            self.s6.set_ydata(np.append(self.s6.get_ydata(),s[5]))
            self.s7.set_ydata(np.append(self.s7.get_ydata(),s[6]))
            self.s8.set_ydata(np.append(self.s8.get_ydata(),s[7]))
            self.s9.set_ydata(np.append(self.s9.get_ydata(),s[8]))
            self.t1.set_ydata(np.append(self.t1.get_ydata(),true[0]))
            self.t2.set_ydata(np.append(self.t2.get_ydata(),true[1]))
            self.t3.set_ydata(np.append(self.t3.get_ydata(),true[2]))
            self.t4.set_ydata(np.append(self.t4.get_ydata(),true[3]))
            self.t5.set_ydata(np.append(self.t5.get_ydata(),true[4]))
            self.t6.set_ydata(np.append(self.t6.get_ydata(),true[5]))
            self.t7.set_ydata(np.append(self.t7.get_ydata(),true[6]))
            self.t8.set_ydata(np.append(self.t8.get_ydata(),true[7]))
            self.t9.set_ydata(np.append(self.t9.get_ydata(),true[8]))

            self.sax1.relim()
            self.sax2.relim()
            self.sax3.relim()
            self.sax4.relim()
            self.sax5.relim()
            self.sax6.relim()
            self.sax7.relim()
            self.sax8.relim()
            self.sax9.relim()
            self.sax1.autoscale_view()
            self.sax2.autoscale_view()
            self.sax3.autoscale_view()
            self.sax4.autoscale_view()
            self.sax5.autoscale_view()
            self.sax6.autoscale_view()
            self.sax7.autoscale_view()
            self.sax8.autoscale_view()
            self.sax9.autoscale_view()
            
            self.sensor_fig.tight_layout()
            self.sensor_fig.canvas.show()

    def ekfplot2(self,s,true,count):
        if self.sensor_first == True:
            # App Thingy
            self.app = pg.QtGui.QApplication([])
            self.app.aboutToQuit.connect(self.stop)

            # Plot Window
            self.ekfplotwin = pg.GraphicsWindow(size=(1200,800))
            self.ekfplotwin.setWindowTitle('Estimated States')
            self.ekfplotwin.setInteractive(True)

            # Plots and subplots
            self.p1 = []
            labels = [u'\u03C6 (rad)',u'\u03B8 (rad)','North (m)','East (m)',
                      'Ground Speed (m/s)',u'\u03C7 (rad)',u'\u03C8 (rad)','','']
            
                      
            for i in xrange(3):
                for j in xrange(3):
                    if i*3 + j < 7:
                        self.p1.append(self.ekfplotwin.addPlot(i+1,j+1))
                        self.p1[i*3 + j].setLabel('bottom','Time',units='s')
                        self.p1[i*3 + j].setLabel('left',labels[i*3 + j])
                        self.p1[i*3 + j].addLegend()
                    
            # Curves to put on plots
            self.estcurves = []
            self.truecurves = []

            # Data for pltos
            self.est = [[],[],[],[],[],[],[],[],[]]
            self.true = [[],[],[],[],[],[],[],[],[]]

            # Put data on curves on correct plots
            for i in xrange(len(s)):
                pen = pg.mkPen(color=(0,255,0),style=pg.QtCore.Qt.DashLine)
                self.estcurves.append(self.p1[i].plot(pen=pen,name='est',))
                pen = pg.mkPen(color=(0,0,255))
                self.truecurves.append(self.p1[i].plot(pen=pen,name='true'))
                self.est[i].append(s[i])
                self.true[i].append(true[i])
            self.sensor_first = False
        else:
            # Update data in plots
            time = np.linspace(0,(len(self.est[0]) + 1)*self.dt,len(self.est[0]) + 1)
            for i in xrange(len(s)):
                self.est[i].append(s[i])
                self.true[i].append(true[i])
                self.estcurves[i].setData(time,self.est[i])
                self.truecurves[i].setData(time,self.true[i])

            # Make it update
            self.app.processEvents()
            
    def stop(self):
        sys.exit()
