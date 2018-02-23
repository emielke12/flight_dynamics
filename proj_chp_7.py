from aircraftClasses.sensors import *
import sys

def ctrl_c(plane,x0,wind):
    # Solution empty vector
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]
    alpha = []
    beta = []
    va = []
    step = []
    counter = 1
    
    try:
        while True:
            # Calculate Wind Values
            # wind[3],wind[4],wind[5] = plane.calc_dryden_gust(plane.dt)

            # Just keep from sideslipping
            plane.sideslip_hold(x0,plane.dt)

            # Calculate Force, Airspeed, alpha, beta
            fx,fy,fz,l,m,n,va_,alpha_,beta_,wn,we,wd = plane.force_calc(x0,plane.deltas,wind)

            # Run ODE to get next step
            sol = odeint(plane.eom,x0,[0.0,plane.dt],args=(fx,fy,fz,l,m,n))

            # Get Sensor Measurements
            if plane.t_sim % plane.Ts < plane.dt:
                sensors = plane.sensor_readings(sol[-1],fx,fy,fz)
            if plane.t_sim % plane.Ts_gps < plane.dt:
                gps = plane.gps_readings(sol[-1],wn,we,wd)

            # Plot Sensors
            plane.plot_sensors_real_time(gps+sensors,counter)

            # Put into total solution matrix for plotting
            for j in xrange(len(x0)):
                # Get height, not z position
                if j == 2:
                    sols[j].append(-sol[-1][j])
                else:
                    sols[j].append(sol[-1][j])

                # Plot chi instead of psi
                if j == 8:
                    sols[j][-1] = plane.chi

            # Alpha Beta and Va Vectors for plotting
            alpha.append(alpha_)
            beta.append(beta_)
            va.append(va_)

            # Plot
            # plane.draw_update([sol[-1][6],sol[-1][7],sol[-1][8]],[sol[-1][0],sol[-1][1],sol[-1][2]])
            # plane.plot_states_real_time(np.transpose(sols))

            # Update initial state for ODE
            x0 = sol[-1,:]

            # Update Time
            plane.t_sim += plane.dt
            counter += 1
            
    except KeyboardInterrupt:
        plt.close('all')
        plane.plot_all_post(np.transpose(sols),alpha,beta,va)
        plt.show()

if __name__ == "__main__":    
    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = 0#-300.0
    u0 = 30.0
    v0 = 0
    w0 = 0.0
    phi0 = 0
    th0 = 0
    psi0 = 0
    p0 = 0.0
    q0 = 0
    r0 = 0

    # Initial State Vector
    x0 = [pn0,pe0,pd0,u0,v0,w0,phi0,th0,psi0,p0,q0,r0]
    
    # Initial Wind Vector
    wind = [0.0,0.0,0.0,0.0,0.0,0.0]

    # Desired Trim
    trim = [x0[3],0.0 * np.pi/180.0,5e2] # Va, gamma, R

    # Instantiate Class
    plane = Sensors(x0,trim)

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(plane,x0,wind)
