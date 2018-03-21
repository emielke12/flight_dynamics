from aircraftClasses.kalmanFilter import *
from aircraftClasses.bottleDrop import *
import sys

def ctrl_c(plane,x0,wind):
    # Solution empty vector
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]
    alpha,beta,va = [],[],[]
    cmd = [[],[],[]]
    counter = 1
    plane.switch = 'air'
    
    try:
        while True:
            # Calculate Wind Values
            # wind[3],wind[4],wind[5] = plane.calc_dryden_gust(plane.dt)

            # Just keep from sideslipping
            # plane.sideslip_hold(x0,plane.dt)
            # plane.sideslip_hold(plane.x_hat,plane.dt)

            ######## Bottle Drop Stuff #######
            plane.wind_estimate(plane.x_hat)
            ######## Bottle Drop Stuff #######


            # Maneuver to Run
            if plane.t_sim >= 40.0:
                plane.altitude_hold(plane.x_hat,105.0,plane.dt)
                plane.airspeed_throttle(plane.x_hat,30.0,plane.dt)
                plane.course_hold(plane.x_hat, 25.0 * np.pi / 180.0, plane.dt,plane.chi_hat)
                cmd[0].append(105.0)
                cmd[1].append(30.0)
                cmd[2].append(25 * np.pi / 180.0)

            else:
                plane.altitude_hold(plane.x_hat,100.0,plane.dt)
                plane.airspeed_throttle(plane.x_hat,30.0,plane.dt)
                plane.course_hold(plane.x_hat, 0.0 * np.pi / 180.0, plane.dt,plane.chi_hat)
                cmd[0].append(100.0)
                cmd[1].append(30.0)
                cmd[2].append(0 * np.pi / 180.0)

                ######## Bottle Drop Stuff #######
                plane.calc_drop_location(-plane.x_hat[2],[100,0],[30.0,0,0])
                ######## Bottle Drop Stuff #######
                
            ######## Bottle Drop Stuff #######
            if plane.in_drop() == True:
                plane.release_triggered(x0,wind)
            ######## Bottle Drop Stuff #######
            
            # Calculate Force, Airspeed, alpha, beta
            fx,fy,fz,l,m,n,va_,alpha_,beta_,wn,we,wd = plane.force_calc(x0,plane.deltas,wind)

            # Run ODE to get next step
            sol = odeint(plane.eom,x0,[0.0,plane.dt],args=(fx,fy,fz,l,m,n))

            # Run EKF
            plane.ekf(sol[-1],[fx,fy,fz],[wn,we,wd],counter,plot = False)
            plt.close('all')

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

            # # Plot
            # plane.draw_update([sol[-1][6],sol[-1][7],sol[-1][8]],[sol[-1][0],sol[-1][1],sol[-1][2]])

            # plot things
            alpha.append(alpha_)
            beta.append(beta_)
            va.append(va_)
            
            # Update initial state for ODE
            x0 = sol[-1,:]

            # Update Time
            plane.t_sim += plane.dt
            counter += 1

            if plane.t_sim >= 40.0:
                plt.close('all')
                plane.plot_all_post(np.transpose(sols),alpha,beta,va,cmd)
                plt.show()
                break
            
    except KeyboardInterrupt:
        plt.close('all')
        sys.exit()

if __name__ == "__main__":    
    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = -100.0#-300.0
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
    wind = [3.0,0.0,0.0,0.0,0.0,0.0]

    # Desired Trim
    trim = [x0[3],0.0 * np.pi/180.0,5e60] # Va, gamma, R

    # Instantiate Class
    plane = bottleDrop(x0,trim)

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(plane,x0,wind)
