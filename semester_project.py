from aircraftClasses.pathPlanner import *
from aircraftClasses.kalmanFilter import *
from aircraftClasses.bottleDrop import *
import sys

def ctrl_c(plane,x0,wind):
    # Solution empty vector
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]
    alpha,beta,va = [],[],[]
    cmd = [[],[],[]]
    counter = 1
    planned = False
    
    try:
        while True:
            # Calculate Wind Values
            # wind[3],wind[4],wind[5] = plane.calc_dryden_gust(plane.dt)

            ######## Bottle Drop Stuff #######
            bottle.wind_estimate(plane.x_hat)
            ######## Bottle Drop Stuff #######

            p = [plane.pn_hat,plane.pe_hat,plane.h_hat]

            # Maneuver to Run
            if plane.t_sim < 1.0:
                r = [0, 0, 100]
                q = [1, 0, 0]
                c = [0, 0, 100]
                rho = 1000
                path_type = 'straight'
                W = [[0,0,100]]
                plan.waypoint_plot(W,p)
            else:
                R = plan.Rmin
                bottle.calc_drop_location(-p[2],[347,275],[30.0, 0, 0,])
                if planned == False:
                    print 'planning'
                    P = plan.planRRT([p[0],p[1],-p[2]], bottle.pdrop, plane.chi_hat)
                    print 'done planning'
                    planned = True
                    P[1][-1] = bottle.approach_angle
                    W = np.transpose(P)[0:3,:]
                    print P

                ######## Bottle Drop Stuff #######
                if bottle.in_drop() == True:
                    bottle.release_triggered(x0,wind)
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

            # Update initial state for ODE
            x0 = sol[-1,:]

            # Update Time
            plane.t_sim += plane.dt
            counter += 1

            plan.waypoint_plot(W,p)

            if plane.t_sim >= 40.0:
#                 plt.close('all')
#                 plane.plot_all_post(np.transpose(sols),alpha,beta,va,cmd)
#                 plt.show()
                break
            
    except KeyboardInterrupt:
        plt.close('all')
        sys.exit()

if __name__ == "__main__":    
    # Initial Conditions
    va0 = 30.0
    x0 = [0,0,-100,va0,0,0,0,0,0,0,0,0]
    model = 1
    
    # Initial Wind Vector
    wind = [3.0,0.0,0.0,0.0,0.0,0.0]

    # Desired Trim
    trim = [x0[3],0.0 * np.pi/180.0,5e60] # Va, gamma, R

    # Instantiate Class
    plan = pathPlanner(va0,model)
    plan.x[4] = -x0[2]
    plane = kalmanFilter(x0,trim)
    plt.close('all')
    bottle = bottleDrop(x0,trim)

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(plane,x0,wind)
