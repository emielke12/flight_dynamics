from autoPilot import *

class Sensors(autoPilot):
    def __init__(self,x0,trim):
        autoPilot.__init__(self,x0,trim)

        # Rates
        self.Ts = 0.02 # seconds
        self.Ts_gps = 1.0 # seconds
        
        # Sensor Variables
        self.sigma_gyro = 0.13 # deg/s
        self.sigma_accel = 0.0025 # m/s**2
        self.beta_abs_pres = 0.125 # kPa
        self.sigma_abs_pres = 0.01 # kPa
        self.beta_diff_pres = 0.02 # kPa
        self.sigma_diff_pres = 0.002 # kPa
        self.sigma_mag = 0.3 # degrees
        self.beta_mag = 1.0 # degrees
        self.sigma_V = 0.05 # m/s

        # GPS Specific Stuff
        self.sigma_n = 0.21 # meters
        self.sigma_e = 0.21 # meters
        self.sigma_d = 0.40 # meters
        self.k_gps = 1.0/1100.0 # Hz
        self.eta_n = 0.0 # m
        self.eta_e = 0.0 # m
        self.eta_d = 0.0 # m

        self.sensor_first = True

    def sensor_readings(self,x,fx,fy,fz):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

        # Rate Gyros
        y_gyro_x = p + np.random.randn() * self.sigma_gyro
        y_gyro_y = q + np.random.randn() * self.sigma_gyro
        y_gyro_z = r + np.random.randn() * self.sigma_gyro

        # Accelerometers
        y_accel_x = fx / self.mass + self.g * np.sin(th) + np.random.randn() * self.sigma_accel
        y_accel_y = fy / self.mass - self.g * np.cos(th) * np.sin(phi) + np.random.randn() * self.sigma_accel
        y_accel_z = fz / self.mass - self.g * np.cos(th) * np.cos(phi) + np.random.randn() * self.sigma_accel

        # Pressure Sensors
        y_abs_pres = self.rho * self.g * -pd + self.beta_abs_pres + np.random.randn() * self.sigma_abs_pres
        y_diff_pres = self.rho * self.Va**2 / 2 + self.beta_diff_pres + np.random.randn() * self.sigma_diff_pres

        return [y_gyro_x,y_gyro_y,y_gyro_z,y_accel_x,y_accel_y,y_accel_z,y_abs_pres,y_diff_pres]

    def gps_readings(self,x,wn,we,wd):
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        
        # Eta update
        self.eta_n = np.exp(-self.k_gps*self.Ts_gps) * self.eta_n + np.random.randn() * self.sigma_n
        self.eta_e = np.exp(-self.k_gps*self.Ts_gps) * self.eta_e + np.random.randn() * self.sigma_e
        self.eta_d = np.exp(-self.k_gps*self.Ts_gps) * self.eta_d + np.random.randn() * self.sigma_d

        # GPS Position Measurements
        y_gps_n = pn + self.eta_n
        y_gps_e = pe + self.eta_e
        y_gps_d = -pd + self.eta_d

        
        # GPS Velocity Measurements
        VN = self.Va * np.cos(psi) + wn
        VE = self.Va * np.sin(psi) + we
        Vg = np.sqrt(VN**2 + VE**2)
        sigma_chi = self.sigma_V / Vg
        y_gps_Vg = np.sqrt(VN**2 + VE**2) + np.random.randn() * self.sigma_V
        y_gps_chi = np.arctan2(VE,VN) + np.random.randn() * sigma_chi

        return [y_gps_n,y_gps_e,y_gps_d,y_gps_Vg,y_gps_chi]

    def plot_sensors_real_time(self,s,count):
        if self.sensor_first == True:
            self.sensor_fig = plt.figure(figsize=(20,10))
            self.sax1 = self.sensor_fig.add_subplot(531)
            self.sax2 = self.sensor_fig.add_subplot(532)
            self.sax3 = self.sensor_fig.add_subplot(533)
            self.sax4 = self.sensor_fig.add_subplot(534)
            self.sax5 = self.sensor_fig.add_subplot(535)
            self.sax6 = self.sensor_fig.add_subplot(536)
            self.sax7 = self.sensor_fig.add_subplot(537)
            self.sax8 = self.sensor_fig.add_subplot(538)
            self.sax9 = self.sensor_fig.add_subplot(539)
            self.sax10 = self.sensor_fig.add_subplot(5,3,10)
            self.sax11 = self.sensor_fig.add_subplot(5,3,11)
            self.sax12 = self.sensor_fig.add_subplot(5,3,12)
            self.sax13 = self.sensor_fig.add_subplot(5,3,13)

            self.sensor_fig.show()
            self.sensor_first = False

            self.sax1.set_ylabel(r'$p_{n,gps}$ (m)')
            self.sax2.set_ylabel(r'$P_{e,gps}$ (m)')
            self.sax3.set_ylabel(r'$h_{gps}$ (m)}')
            self.sax4.set_ylabel('V_g (m/s)')
            self.sax5.set_ylabel(r'$\chi$ (rad)')
            self.sax6.set_ylabel(r'$gyro_x$ (rad/s)')
            self.sax7.set_ylabel(r'$gyro_y$ (rad/s)')
            self.sax8.set_ylabel(r'$gyro_z$ (rad/s)')
            self.sax9.set_ylabel(r'$a_x$ (rad/s)')
            self.sax10.set_ylabel(r'$a_y$ (rad/s)')
            self.sax11.set_ylabel(r'$a_z$ (rad/s)')
            self.sax12.set_ylabel(r'$P_{abs}$ (kPa)')
            self.sax13.set_ylabel(r'$P_{diff}$ (kPa)')

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
            self.s12, = self.sax12.plot(s[11])
            self.s13, = self.sax13.plot(s[12])

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
            self.s12.set_xdata(t)
            self.s13.set_xdata(t)

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
            self.s12.set_ydata(np.append(self.s12.get_ydata(),s[11]))
            self.s13.set_ydata(np.append(self.s13.get_ydata(),s[12]))

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
            self.sax12.relim()
            self.sax13.relim()
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
            self.sax12.autoscale_view()
            self.sax13.autoscale_view()

            
            self.sensor_fig.tight_layout()
            self.sensor_fig.canvas.draw()


        
