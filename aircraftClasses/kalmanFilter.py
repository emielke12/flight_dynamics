from sensors import *
from copy import deepcopy

class kalmanFilter(Sensors):
    def __init__(self,x0,trim):
        Sensors.__init__(self,x0,trim)
        self.sens = [0.0] * 8
        self.gps = [0.0] * 5
        self.Tout = 0.01
        
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
        self.Vg_hat = 0.0

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
        self.alpha_pn = 0.2
        self.alpha_pe = 0.2
        self.alpha_chi = 0.2
        self.alpha_vg = 0.2

        # Previous Sensor Values
        self.u_att = [0.0,0.0,0.0,0.0]

    def ekf(self,x_att,x_gps,u_att,u_gps):
        # roll and pitch stuff
        x_atthat,P_att = predict(x_att,self.f_att,u_att,self.df_att,self.Pa,self.Qa)
        for i in xrange(len(u_att)):
            if u_att[i] != self.u_att[i]:
                x_atthat,P_att = update(x_att_hat,self.dh_att,u_att,P_att,self.Ra,y[i],self.h_att,i)

        self.Pa = P_att
        self.x_atthat = x_atthat

        # Gps stuff

    def predict(self,xhat,f,u,df,P,Q):
        for i in xrange(10):
            xhat += self.Tout/i * f(xhat,u)
            A = df(xhat,u)
            P += self.Tout/i * (A * P + P * A.T + Q)
        return xhat,P

    def update(self,xhat,dh,u,P,R,y,h,col):
        C = dh(xhat[col],u)
        L = P * C.T * 1 / (R + C * P * C.T)
        P = (np.eye(np.shape(L * C)) - L * C) * P
        xhat += L * (y - h(xhat,u))
        return xhat,P,L

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
        f = [[q * np.cos(phi) * np.tan(th) - r * np.sin(phi) * np.tan(th), (q * np.sin(phi) - r * np.cos(phi)) / (np.cos(th)**2)],
             [-q * np.sin(phi) - r * np.cos(phi), 0]]
        return f

    def dh_att(self,x,u,col):
        phi, th = x
        p, q, r, va = u
        f = [[0, q * va * np.cos(th) + self.g * np.cos(th)],
             [-self.g * np.cos(phi) * np.cos(th), -r * va * np.sin(th) - p * va * np.cos(th) + self.g * np.sin(phi) * np.sin(th)],
             [self.g * np.sin(phi) * np.cos(th), (q * va + self.g * np.cos(phi)) * np.sin(th)]]
        return f[col]
    
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
        self.p_hat = self.lpf(self.p_hat,y_gyro_x,self.alpha_p)
        self.q_hat = self.lpf(self.q_hat,y_gyro_y,self.alpha_q)
        self.r_hat = self.lpf(self.r_hat,y_gyro_z,self.alpha_r)
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
        self.phi_hat = np.arctan(self.ay_hat / self.az_hat)
        self.th_hat = np.arcsin(self.ax_hat / self.g)
        self.pn_hat = self.lpf(self.pn_hat, y_gps_n, self.alpha_pn)
        self.pe_hat = self.lpf(self.pe_hat, y_gps_e, self.alpha_pe)
        self.chi_hat = self.lpf(self.chi_hat, y_gps_chi, self.alpha_chi)
        self.Vg_hat = self.lpf(self.Vg_hat, y_gps_Vg, self.alpha_vg)

        # Plot
        s = [self.pn_hat,self.pe_hat,self.h_hat,self.Vg_hat,self.chi_hat,
             self.p_hat,self.q_hat,self.r_hat,self.Va_hat,self.phi_hat,self.th_hat]
        true = [pn,pe,-pd,self.Vg,self.chi,p,q,r,self.Va,phi,th]
        print s[2] + pd
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
