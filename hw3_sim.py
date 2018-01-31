from aircraftClasses.force_moments import *
import sys

if __name__ == "__main__":
    meq = MAVForces()

    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = 0
    u0 = 10.0
    v0 = 0
    w0 = 0.0
    phi0 = 0
    th0 = 0
    psi0 = 0
    p0 = 0.0
    q0 = 0
    r0 = 0
    meq.Va = np.linalg.norm([u0,v0,w0])

    x0 = [pn0,pe0,pd0,u0,v0,w0,phi0,th0,psi0,p0,q0,r0]

    # Time Steps
    step = 300

    # # Inputs
    fx = np.zeros(step)
    fy = np.zeros(step)
    fz = np.zeros(step)
    l = np.zeros(step)
    m = np.zeros(step)
    n = np.zeros(step)

    # Control Surface Inputs and Wind
    wind = [0.0,0.0,0.0,0.0,0.0,0.0] # steady wind is last 3
    cs = np.zeros((4,step)) # elevator,aileron,rudder,thrust
    cs[0][:] = -10.576 * np.pi/180.0 
    cs[3][:] = 0.3
    cs[2][:] = 0.0

    # Time Vector
    t = np.linspace(0,10,step + 1)

    # Solve ODE
    try:
        sols = np.zeros((len(t)-1,len(x0)))
        for i in xrange(step):
            # Calculate Wind Values
            wind[3],wind[4],wind[5] = meq.calc_dryden_gust(t[1] - t[0])

            # Calculate Force, Airspeed, alpha, beta
            fx,fy,fz,l,m,n,va,alpha,beta,wn,we,wd = meq.force_calc(x0,cs[:,i],wind)
            wind[0] = wn
            wind[1] = we
            wind[2] = wd

            # Run ODE to get next step
            sol = odeint(meq.eom,x0,[t[i],t[i+1]],args=(fx,fy,fz,l,m,n))

            # Plot
            meq.plot_update([sol[-1][6],sol[-1][7],sol[-1][8]],[sol[-1][0],sol[-1][1],sol[-1][2]])

            for j in xrange(len(x0)):
                sols[i,j] = sol[-1][j]

            # Update initial state for ODE
            x0 = sol[-1:][0]
    except KeyboardInterrupt:
        sys.exit()

    meq.plot_states(sols[:i+1,:])
    # plt.close('all')
    # fig = plt.figure()
    # ax1 = fig.add_subplot(261)
    # ax1.plot(sols[:,0])
    # ax2 = fig.add_subplot(262)
    # ax2.plot(sols[:,1])
    # ax3 = fig.add_subplot(263)
    # ax3.plot(sols[:,2])
    # ax4 = fig.add_subplot(264)
    # ax4.plot(sols[:,3])
    # ax5 = fig.add_subplot(265)
    # ax5.plot(sols[:,4])
    # ax6 = fig.add_subplot(266)
    # ax6.plot(sols[:,5])
    # ax7 = fig.add_subplot(267)
    # ax7.plot(sols[:,6])
    # ax8 = fig.add_subplot(268)
    # ax8.plot(sols[:,7])
    # ax9 = fig.add_subplot(269)
    # ax9.plot(sols[:,8])
    # ax10 = fig.add_subplot(2,6,10)
    # ax10.plot(sols[:,9])
    # ax11 = fig.add_subplot(2,6,11)
    # ax11.plot(sols[:,10])
    # ax12 = fig.add_subplot(2,6,12)
    # ax12.plot(sols[:,11])
    # plt.show()
    # pos = [sols[:,0],sols[:,1],sols[:,2]]
    # eul = [sols[:,6],sols[:,7],sols[:,8]]
    # meq.input_goal(eul,pos)
