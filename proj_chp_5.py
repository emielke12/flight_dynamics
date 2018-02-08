from aircraftClasses.mavTrim import *
import sys

def ctrl_c(meq,x0):
    # Initial Wind Values
    wind = [0.0,0.0,0.0,0.0,0.0,0.0]

    # Desired Trim
    trim = [x0[3],0.0 * np.pi/180.0,5e60] # Va, gamma, R

    # Solution empty vector
    sols = [[],[],[],[],[],[],[],[],[],[],[],[]]

    # Initial guess for alpha,beta,phi
    star = [0.0,0.0,0.0]

    # Optimize Deltas
    meq.call_opt(star,x0,trim,wind)
    x0 = meq.update_init_conds
    # print 'Initial Conditions: ',x0,'\n'
    # print 'Elevator, Aileron, Rudder, Thrust: ',meq.deltas

    # Get Transfer Functions
    meq.get_transfer_functions(x0)
    # print meq.T_phi_deltaa
    
    # Get SS Matrices
    meq.get_numerical_ss(x0)
    # print 'A_lon:\t',meq.A_lon,'\n'
    # print 'B_lon:\t',meq.B_lon,'\n'
    # print 'Eig A_lon:\t',meq.eig_lon,'\n'
    # print 'Eig A_lat:\t',meq.eig_lat,'\n'

    # raw_input()
    
    # print meq.A_lat
    # print meq.B_lat
    # raw_input()


    try:
        counter = 0
        use = 'dutch roll'
        while True:
            # Calculate Wind Values
            # wind[3],wind[4],wind[5] = meq.calc_dryden_gust(0.03)

            # Input Impulse or Doublet
            if use == 'phugoid':
                if counter == 100:
                    imp = [-5.0 * np.pi / 180.0, 0, 0, 0]
                    meq.deltas = np.add(meq.deltas,imp).tolist()
                elif counter == 116:
                    imp = [5.0 * np.pi / 180.0, 0, 0, 0]
                    meq.deltas = np.add(meq.deltas,imp).tolist()
            elif use == 'dutch roll':
                if counter == 100:
                    doub = [0, 0, 25 * np.pi / 180.0, 0]
                    meq.deltas = np.add(meq.deltas,doub).tolist()    
                elif counter == 108:
                    doub = [0, 0, -50 * np.pi / 180.0, 0]
                    meq.deltas = np.add(meq.deltas,doub).tolist()
                elif counter == 116:
                    doub = [0, 0, 25 * np.pi / 180.0, 0]
                    meq.deltas = np.add(meq.deltas,doub).tolist()
            elif use == 'spiral':
                if counter == 100:
                    doub = [0, 25 * np.pi / 180.0, 0, 0]
                    meq.deltas = np.add(meq.deltas,doub).tolist()    
                elif counter == 108:
                    doub = [0, -50 * np.pi / 180.0, 0, 0]
                    meq.deltas = np.add(meq.deltas,doub).tolist()
                elif counter == 116:
                    doub = [0, 25 * np.pi / 180.0, 0, 0]
                    meq.deltas = np.add(meq.deltas,doub).tolist()
            else:
                pass


            # Calculate Force, Airspeed, alpha, beta
            fx,fy,fz,l,m,n,va,alpha,beta,wn,we,wd = meq.force_calc(x0,meq.deltas,wind)

            # Run ODE to get next step
            sol = odeint(meq.eom,x0,[0.0,0.03],args=(fx,fy,fz,l,m,n))
            # sol[:,6] = (sol[:,6] + np.pi) % (2 * np.pi) - np.pi
            # sol[:,7] = (sol[:,7] + np.pi) % (2 * np.pi) - np.pi
            # sol[:,8] = (sol[:,8] + np.pi) % (2 * np.pi) - np.pi

            # Put into total solution matrix
            for j in xrange(len(x0)):
                # Get height, not z position
                if j == 2:
                    sols[j].append(-sol[-1][j])
                else:
                    sols[j].append(sol[-1][j])

            # Plot
            meq.draw_update([sol[-1][6],sol[-1][7],sol[-1][8]],[sol[-1][0],sol[-1][1],sol[-1][2]])
            # meq.plot_states_real_time(np.transpose(sols))

            # Update initial state for ODE
            x0 = sol[-1,:]

            # Update counter
            counter += 1
            
    except KeyboardInterrupt:
        plt.close('all')
        meq.plot_states_post(np.transpose(sols))

if __name__ == "__main__":
    # Instantiate Class
    meq = Trim()
    
    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = 0#-300.0
    u0 = 17.0
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
