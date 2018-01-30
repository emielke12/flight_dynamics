from aircraftClasses.force_moments import *

if __name__ == "__main__":
    meq = MAVForces()

    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = 0
    u0 = 1.0
    v0 = 0
    w0 = 1.0
    phi0 = 0
    th0 = 0
    psi0 = 0
    p0 = -1.0
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
    wind = [0.0,0.0,0.0,0.1,0.1,-0.1]
    cs = np.zeros((4,step)) # elevator,aileron,rudder,thrust
    cs[0][:] = 0.0
    cs[3][0:150] = 0.5
    cs[2][:] = 0.1

    # Time Vector
    t = np.linspace(0,10,step + 1)

    # Solve ODE
    sols = np.zeros((len(t)-1,len(x0)))
    for i in xrange(step):
        wind[0],wind[1],wind[2] = meq.calc_dryden_gust(t[1] - t[0])
        fx,fy,fz,l,m,n,va,alpha,beta,wn,we,wd = meq.force_calc(x0,cs[:,i],wind)
        sol = odeint(meq.eom,x0,[t[i],t[i+1]],args=(fx,fy,fz,l,m,n))
        for j in xrange(len(x0)):
            sols[i,j] = sol[-1][j]
        x0 = sol[-1:][0]        


    pos = [sols[:,0],sols[:,1],sols[:,2]]
    eul = [sols[:,6],sols[:,7],sols[:,8]]
    meq.input_goal(eul,pos)
