from aircraftClasses.pathPlanner import *
from aircraftClasses.bottleDrop import *
import sys

def ctrl_c(plane,plane2,bottle,x0,wind,model):
    # Kalman Filter AutoPilot emtpy vectors
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]

    # command for kalman autopilot
    cmd2 = [[],[],[]]
    counter = 1

    # Command for simplified 1
    cmd = [0,0,0,100,plane.va0,0] # chidot,chi,hdot,h,va,phi_c

    # Use chapter 9 or kalman filter
    kalman = True
    path_type = 'straight'
    switch_time = 2.0

    try:
        while True:
            # Position and commanded Velocity
            p_9 = [plane.x[0],plane.x[1],plane.x[4]]
            p_est = [plane2.pn_hat,plane2.pe_hat,plane2.h_hat]
            cmd[4] = 30.0 # va

            # Estimate wind
            bottle.wind_estimate(plane2.x_hat,plane2.chi_hat,plane2.psi_hat,wind)

            # Maneuver to Run
            if plane.t_sim < switch_time:
                r = [0, 0, 100]
                q = [1, 0, 0]
                c = [0,0,100]
                rho = 1000
                path_type = 'straight'
                W = [[0,0,100]]

            elif plane.t_sim >= switch_time and plane.t_sim < switch_time + 0.01:
                R = plane.Rmin
                target = [600,400]
#                 bottle.calc_drop_location(-p_9[2],[670,170],[30,0,0])
                # bottle.calc_drop_location(-cmd[3],target,[30,0,0],False)
                bottle.calc_drop_location(-cmd[3],target,[0,-30,0],True)
                bottle.pdrop[2] = -100.0
#                 P = plane.planRRT([p_9[0], p_9[1], -p_9[2]],bottle.pdrop,plane.x[2])
                # P = plane.planRRT([p_est[0], p_est[1], -p_est[2]],bottle.pdrop,plane2.chi_hat)
                # P[1][-1] = bottle.approach_angle
                P = bottle.plan_path([p_est[0], p_est[1], -p_est[2]],bottle.pdrop,plane2.chi_hat)

                W = np.transpose(np.transpose(P)[0:3,:])
#                 path_type,r,q,c,rho,plane.lamb = plane.algorithm_8(P,p_9,R)
#                 plane.plot_arrow(plane.x[2],p_9)
                path_type,r,q,c,rho,plane.lamb = plane.algorithm_8(P,p_est,R)
                plane.plot_arrow(plane2.chi_hat,p_est)

            elif plane.t_sim >= switch_time + 0.01:
                R = plane.Rmin
#                 path_type,r,q,c,rho,plane.lamb = plane.algorithm_8(P,p_9,R)
#                 plane.plot_arrow(plane.x[2],p_9)
                path_type,r,q,c,rho,plane.lamb = plane.algorithm_8(P,p_est,R)
                plane.plot_arrow(plane2.chi_hat,p_est)

                if bottle.in_drop(p_est) == True:
                    bottle.release_triggered(x0,wind,target)

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
            # print plane2.chi_hat - cmd2[2][-1]

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
            # plane.plot(commands,states,autostates,labels)
            if kalman == True:
                plane.waypoint_plot(W,p_est)
            else:
                plane.waypoint_plot(W,p_9)
            
            # Update initial state
            x0 = sol[-1,:]

    except KeyboardInterrupt:
        sys.exit()

if __name__ == "__main__":    
    # Initial Conditions
    va0 = 30.0
    x0 = [0,0,-100,va0,0,0,0,0,0,0,0,0]
    model = 1
    
    # Initial Wind Vector
    # wind = [3,-3,0,0,0,0]
    # wind = [0,1,0,0,0,0]
    wind = [0,0,0,0,0,0]

    # Trim
    trim = [x0[3],0,5e60]
    
    # Instantiate Class
    plane = pathPlanner(va0,model)
    plane.x[4] = -x0[2]
    plane2 = kalmanFilter(x0,trim)
    plt.close('all')
    bottle = bottleDrop(x0,trim)
    plt.close('all')

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(plane,plane2,bottle,x0,wind,model)
