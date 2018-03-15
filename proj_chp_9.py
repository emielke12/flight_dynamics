from aircraftClasses.guidance import *
import sys

def ctrl_c(plane,plane2,x0,wind,model):
    # Kalman Filter AutoPilot emtpy vectors
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]

    # command for kalman autopilot
    cmd2 = [[],[],[]]
    counter = 1

    # Command for simplified 1
    cmd = [0,0,0,0,plane.va0,0] # chidot,chi,hdot,h,va,phi_c

    try:
        while True:
            # Maneuver to Run
            if plane.t_sim >= 1.0:
                cmd[3] = 10.0 # height
                cmd[1] = 20 * np.pi / 180.0 # chi
                cmd[4] = 35.0 # va
                cmd[5] = 30 * np.pi / 180.0 # phi

                plane2.altitude_hold(x0,cmd[3],plane2.dt)
                plane2.airspeed_throttle(x0,cmd[4],plane2.dt)
                cmd2[0].append(cmd[3])
                cmd2[1].append(cmd[4])
                
                if model == 1:
                    plane2.course_hold(x0,cmd[1],plane2.dt)
                    cmd2[2].append(cmd[1])
                elif model == 2:
                    plane2.roll_attitude(x0,cmd[5],plane2.dt)
                    cmd2[2].append(cmd[5])
                
            # Run ODE
            plane.simulate_states(wind[0:3],cmd,model)
            fx,fy,fz,l,m,n,va_,alpha_,beta_,wn,we,wd = plane2.force_calc(x0,plane2.deltas,wind)
            sol = odeint(plane2.eom,x0,[0.0,plane2.dt],args=(fx,fy,fz,l,m,n))
            plane2.ekf(sol[-1],[fx,fy,fz],[wn,we,wd],counter,plot = True)

            # Update Time
            plane.t_sim += plane.dt
            plane2.t_sim += plane2.dt
            counter += 1

            # Plot
            if model == 1:
                commands = [cmd[1],cmd[3],cmd[4]]
                states = [plane.x[2],plane.x[4],plane.x[6]]
                autostates = [sol[-1][8],-sol[-1][2],va_]
                labels = [u'\u03C7 (rad)','h (m)','Va (m/s)'] # chi, height, airspeed
                plane.plot(commands,states,autostates,labels)
            elif model == 2:
                commands = [cmd[3],cmd[4],cmd[5]]
                states = [plane.x[3],plane.x[5],plane.x[6]]
                autostates = [-sol[-1][2],va_,sol[-1][6]]
                labels = ['h (m)','Va (m/s)',u'\u03C6 (rad)'] # Height, airspeed, phi
                plane.plot(commands,states,autostates,labels)

            # Update initial state
            x0 = sol[-1,:]
            
    except KeyboardInterrupt:
        sys.exit()

if __name__ == "__main__":    
    # Initial Conditions
    va0 = 30.0
    x0 = [0,0,0,va0,0,0,0,0,0,0,0,0]
    model = 1
    
    # Initial Wind Vector
    wind = [0.0,0.0,0.0]
    wind2 = [0,0,0,0,0,0]

    # Trim
    trim = [x0[3],0,5e60]
    
    # Instantiate Class
    plane = guidanceModel(va0,model)
    plane2 = kalmanFilter(x0,trim)

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(plane,plane2,x0,wind2,model)
