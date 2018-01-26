from aircraftClasses.force_moments import *

if __name__ == "__main__":
    meq = MAVForces()

    # Initial Conditions
    pn0 = 0
    pe0 = 0
    pd0 = 0
    u0 = 0
    v0 = 0
    w0 = 0
    phi0 = 0
    th0 = 0
    psi0 = 0
    p0 = 0
    q0 = 0
    r0 = 0

    x0 = [pn0,pe0,pd0,u0,v0,w0,phi0,th0,psi0,p0,q0,r0]

    # Time Steps
    step = 421

    # Inputs
    fx = np.zeros(step)
    fy = np.zeros(step)
    fz = np.zeros(step)
    l = np.zeros(step)
    m = np.zeros(step)
    n = np.zeros(step)
    fx[20:40] = 40.0
    fx[40:80] = -40.0
    fx[80:100] = 40.0
    fy[100:120] = 40.0
    fy[120:160] = -40.0
    fy[160:180] = 40.0
    fz[180:200] = 40.0
    fz[200:240] = -40.0
    fz[240:260] = 40.0
    m[260:280] = 10.0
    m[280:320] = -10.0
    m[320:340] = 10.0
    n[340:360] = 10.0
    n[360:400] = -10.0
    n[400:420] = 10.0

    # Control Surface Inputs and Wind
    

    # Time Vector
    t = np.linspace(0,10,step)

    # Solve ODE
    sols = np.zeros((len(t),len(x0)))
    for i in xrange(step):
        sol = odeint(meq.eom,x0,[t[i-1],t[i]],args=(fx[i],fy[i],fz[i],l[i],m[i],n[i]))
        for j in xrange(len(x0)):
            sols[i,j] = sol[-1][j]
        x0 = sol[-1:][0]        


    pos = [sols[:,0],sols[:,1],sols[:,2]]
    eul = [sols[:,6],sols[:,7],sols[:,8]]
    meq.input_goal(eul,pos)
