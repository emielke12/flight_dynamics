from aircraftClasses.autoPilot import *
import sys
from copy import deepcopy

def ctrl_c(meq,x0,wind):
    # Solution empty vector
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]
    alpha = []
    beta = []
    va = []
    step = []

    th_d = meq.update_init_conds[7]

    # # Desired Roll Angle
    # phi_d = 10.0 * np.pi/180.0

    # # Desired Pitch Angle
    # th_d = 10.0 * np.pi / 180.0
    
    # # Desired Course Angle
    # chi_d = 5.0 * np.pi / 180.0

    # # Desired Airspeed
    # va_d = 40.0

    # # Desired Altitude
    # h_d = 20.0
    
    counter = 0

    try:
        while True:
            # Calculate Wind Values
            # wind[3],wind[4],wind[5] = meq.calc_dryden_gust(meq.dt)

            if counter > 20:                
                # Inner Loops
                meq.sideslip_hold(x0,meq.dt)
                # meq.roll_attitude(x0,phi_d,meq.dt)
                # meq.airspeed_pitch(x0,va_d,meq.dt)
                
                # # Course Step
                # meq.course_hold(x0,chi_d,meq.dt)
                # step.append(chi_d)
                
                # # Airspeed Step
                # meq.airspeed_throttle(x0,va_d,meq.dt)
                # meq.altitude_hold(x0,0,meq.dt)
                # step.append(va_d)

                # # Altitude Step
                # meq.altitude_hold(x0,h_d,meq.dt)
                # step.append(h_d)

                # # AutoPilot
                # meq.mode(x0)
            else:
                step.append(0)
                # step.append(meq.Va_0)

            # Calculate Force, Airspeed, alpha, beta
            fx,fy,fz,l,m,n,va_,alpha_,beta_,wn,we,wd = meq.force_calc(x0,meq.deltas,wind)

            # Run ODE to get next step
            sol = odeint(meq.eom,x0,[0.0,meq.dt],args=(fx,fy,fz,l,m,n))

            # Put into total solution matrix for plotting
            for j in xrange(len(x0)):
                # Get height, not z position
                if j == 2:
                    sols[j].append(-sol[-1][j])
                else:
                    sols[j].append(sol[-1][j])

                # Plot chi instead of psi
                if j == 8:
                    sols[j][-1] = meq.chi

            # Alpha Beta and Va Vectors for plotting
            alpha.append(alpha_)
            beta.append(beta_)
            va.append(va_)

            # Plot
            meq.draw_update([sol[-1][6],sol[-1][7],sol[-1][8]],[sol[-1][0],sol[-1][1],sol[-1][2]])
            # meq.plot_states_real_time(np.transpose(sols))

            # Update initial state for ODE
            x0 = sol[-1,:]

            # Update counter
            counter += 1
            
    except KeyboardInterrupt:
        plt.close('all')
        meq.plot_all_post(np.transpose(sols),alpha,beta,va)
        time = np.multiply(np.linspace(0,len(step),len(step)),meq.dt)
        meq.ax3.plot(time,step,color='k',linestyle='--') # Alt
        # meq.ax7.plot(time,step,color='k',linestyle='--') # Phi
        # meq.ax8.plot(time,step,color='k',linestyle='--') # Theta
        # meq.ax9.plot(time,step,color='k',linestyle='--') # Chi
        # meq.ax14.plot(time,step,color='k',linestyle='--') # Beta
        # meq.ax15.plot(time,step,color='k',linestyle='--') # Air
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
    trim = [x0[3],0.0 * np.pi/180.0,5e60] # Va, gamma, R

    # Instantiate Class
    meq = autoPilot(x0,trim)

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(meq,x0,wind)
