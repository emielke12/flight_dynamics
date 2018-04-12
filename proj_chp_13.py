from aircraftClasses.cameras import *
import sys,time

def ctrl_c(plane,plane2,x0,wind,model):
    # Kalman Filter AutoPilot emtpy vectors
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]

    # command for kalman autopilot
    cmd2 = [[],[],[]]
    counter = 1

    # Command for simplified 1
    cmd = [0,0,0,100,plane.va0,0] # chidot,chi,hdot,h,va,phi_c

    # Use chapter 9 or kalman filter
    kalman = False
    path_type = 'straight'

    try:
        while True:
            # Position and commanded Velocity
            p_9 = [plane.x[0],plane.x[1],plane.x[4]]
            p_est = [plane2.pn_hat,plane2.pe_hat,plane2.h_hat]
            # p_est = [x0[0],x0[1],-x0[2]]
            cmd[4] = 30.0 # va
            
            # Maneuver to Run
            r = [0, 0, 100]
            q = [1, 0, 0]
            c = [0,0,100]
            rho = 1000
            path_type = 'straight'
            W = [[0,0,100]]

            # Target Tracking
            p_obj = [10,0,0]
            p_mav = x0[0:3]
            phi,th,psi = x0[6:9]
            plane.geolocation(p_obj,p_mav,phi,th,psi)
            time.sleep(0.1)

            # Determine Commands for autopilot
            if kalman == True:
                if path_type == 'straight':
                    cmd[3], cmd[1] = plane.follow_straight_line(r,q,p_est,plane2.chi_hat)
                    # cmd[3], cmd[1] = plane.follow_straight_line(r,q,p_est,plane2.chi)
                elif path_type == 'orbit':
                    cmd[3], cmd[1] = plane.follow_orbit(c,rho,p_est,plane2.chi_hat)
                    # cmd[3], cmd[1] = plane.follow_orbit(c,rho,p_est,plane2.chi)
            else:
                if path_type == 'straight':
                    cmd[3], cmd[1] = plane.follow_straight_line(r,q,p_9,plane.x[2])
                elif path_type == 'orbit':
                    cmd[3], cmd[1] = plane.follow_orbit(c,rho,p_9,plane.x[2])

            # Keep track of commanded course angle, height, and airspeed
            cmd2[0].append(cmd[3])
            cmd2[1].append(cmd[4])
            cmd2[2].append(cmd[1])

            # Autopilot Kalman Filter states
            plane2.altitude_hold(plane2.x_hat,cmd[3],plane2.dt)
            plane2.airspeed_throttle(plane2.x_hat,cmd[4],plane2.dt)
            plane2.course_hold(plane2.x_hat,cmd[1],plane2.dt,plane2.chi_hat)

            # Run ODE
            plane.simulate_states(wind[0:3],cmd,model)
            fx,fy,fz,l,m,n,va_,alpha_,beta_,wn,we,wd = plane2.force_calc(x0,plane2.deltas,wind)
            sol = odeint(plane2.eom,x0,[0.0,plane2.dt],args=(fx,fy,fz,l,m,n))
            plane2.ekf(sol[-1],[fx,fy,fz],[wn,we,wd],counter,plot = False)

            # Update Time
            plane.t_sim += plane.dt
            plane2.t_sim += plane2.dt
            counter += 1

            # Plot
            commands = [cmd[1],cmd[3],cmd[4]]
            states = [plane.x[2],plane.x[4],plane.x[6]]
            autostates = [plane2.chi_hat,plane2.h_hat,plane2.Va_hat]
            labels = [u'\u03C7 (rad)','h (m)','Va (m/s)'] # chi, height, airspeed


#             if kalman == True:
#                 plane.waypoint_plot(W,p_est)
#             else:
#                 plane.waypoint_plot(W,p_9)
            
            # Update initial state
            x0 = sol[-1,:]

    except KeyboardInterrupt:
        print '\nGeolocation Error:\t',plane.geo_error
        sys.exit()

if __name__ == "__main__":    
    # Initial Conditions
    va0 = 30.0
    x0 = [0,0,-100,va0,0,0,0,0,0,0,0,0]
    model = 1
    
    # Initial Wind Vector
    # wind = [3,-3,0,0,0,0]
    wind = [0,0,0,0,0,0]

    # Trim
    trim = [x0[3],0,5e60]
    
    # Instantiate Class
    plane = Cameras(va0,model)
    plane.x[4] = -x0[2]
    plane2 = kalmanFilter(x0,trim)
    plt.close('all')

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(plane,plane2,x0,wind,model)
