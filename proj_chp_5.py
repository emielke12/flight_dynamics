from aircraftClasses.mavTrim import *
import sys

def ctrl_c(meq,x0):
    # Initial Wind Values
    wind = [0.0,0.0,0.0,0.0,0.0,0.0]

    # Desired Trim
    trim = [x0[3],-5.0 * np.pi/180.0,1000000000000000] # Va, gamma, R

    # Solution empty vector
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]

    # Initial guess for alpha,beta,phi
    star = [0.0,0.0,0.0]

    # Optimize Deltas
    meq.call_opt(star,x0,trim,wind)
    x0 = meq.update_init_conds

    try:
        while True:
            # Calculate Wind Values
            # wind[3],wind[4],wind[5] = meq.calc_dryden_gust(0.03)

            # Calculate Force, Airspeed, alpha, beta
            fx,fy,fz,l,m,n,va,alpha,beta,wn,we,wd = meq.force_calc(x0,meq.deltas,wind)

            # Run ODE to get next step
            sol = odeint(meq.eom,x0,[0.0,0.03],args=(fx,fy,fz,l,m,n))
            sol[:,6] = (sol[:,6] + np.pi) % (2 * np.pi) - np.pi
            sol[:,7] = (sol[:,7] + np.pi) % (2 * np.pi) - np.pi
            sol[:,8] = (sol[:,8] + np.pi) % (2 * np.pi) - np.pi

            # Put into total solution matrix
            for j in xrange(len(x0)):
                if j == 2:
                    sols[j].append(-sol[-1][j])
                else:
                    sols[j].append(sol[-1][j])

            # Plot
            meq.draw_update([sol[-1][6],sol[-1][7],sol[-1][8]],[sol[-1][0],sol[-1][1],sol[-1][2]])
            # meq.plot_states_real_time(np.transpose(sols))

            # Update initial state for ODE
            x0 = sol[-1:][0]
            
    except KeyboardInterrupt:
        plt.close('all')
        meq.plot_states_post(np.transpose(sols))

if __name__ == "__main__":
    # Instantiate Class
    meq = Trim()
    
    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = 0
    u0 = 20.0
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

    #------------------------ RUNS UNTIL CTRL-C---------------------------------------#
    ctrl_c(meq,x0)

    #------------------------ RUNS UNTIL STEPS---------------------------------------#
    # # Time Steps
    # step = 300

    # # Control Surface Inputs and Wind
    # wind = [0.0,0.0,0.0,0.0,0.0,0.0] # [0-2] Steady Wind [3-5] Gusts
    # cs = np.zeros((4,step)) # elevator,aileron,rudder,thrust
    # cs[0][:] = -0.567 * np.pi/180.0 
    # cs[3][:] = 0.3
    # cs[2][:] = 0.0
    
    # # Time Vector
    # t = np.linspace(0,10,step + 1)

    # # Solve ODE
    # sols = np.zeros((len(t)-1,len(x0)))

    # try:
    #     for i in xrange(step):
    #         # Calculate Wind Values
    #         wind[3],wind[4],wind[5] = meq.calc_dryden_gust(t[1] - t[0])

    #         # Calculate Force, Airspeed, alpha, beta
    #         fx,fy,fz,l,m,n,va,alpha,beta,wn,we,wd = meq.force_calc(x0,cs[:,i],wind)

    #         # Run ODE to get next step
    #         sol = odeint(meq.eom,x0,[t[i],t[i+1]],args=(fx,fy,fz,l,m,n))

    #         # Put into total solution matrix
    #         for j in xrange(len(x0)):
    #             sols[i,j] = sol[-1][j]

    #         # Plot
    #         # meq.draw_update([sol[-1][6],sol[-1][7],sol[-1][8]],[sol[-1][0],sol[-1][1],sol[-1][2]])
    #         meq.plot_states_real_time(sols[:i,:])

    #         # Update initial state for ODE
    #         x0 = sol[-1:][0]
            
    # except KeyboardInterrupt:
    #     sys.exit()

    # plt.close('all')
    # meq.plot_states_post(sols)



    
