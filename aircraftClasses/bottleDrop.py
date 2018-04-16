from kalmanFilter import *
import pyqtgraph as pg

class bottleDrop(kalmanFilter):
    def __init__(self,x0,trim):
        kalmanFilter.__init__(self,x0,trim)
        self.wind_array = []
        self.wind_avg = [0.0,0.0]
        self.wind_psi = 0.0

        # If object is sphere and low speed
        self.radius = 0.25 # m
        self.mass_object = 0.5 # kg
        self.S_object = np.pi * self.radius ** 2
        self.Cd_object = 0.5

    def wind_estimate(self,x0,chi,psi,wind):
        # Estimate wind from current state or state estimate (if using kalman filter)
        wn = self.Vg_hat * np.cos(chi) - self.Va_hat * np.cos(psi)
        we = self.Vg_hat * np.sin(chi) - self.Va_hat * np.sin(psi)

        if len(self.wind_array) < 50:
            self.wind_array.append([wn,we])
        else:
            self.wind_array = np.roll(self.wind_array,-1,axis=0)
            self.wind_array[-1] = [wn,we]
        self.wind_avg = np.average(self.wind_array,axis=0)
        self.wind_psi = np.arctan2(self.wind_avg[1],self.wind_avg[0]) - np.pi/2
        # print self.wind_psi
        # self.wind_avg = np.add(wind[0:2],np.random.randn(2) * 0.01)
        self.wind_avg = np.array([0.0,0.0])
        # self.wind_psi = np.arctan2(self.wind_avg[1],self.wind_avg[0])
        self.wind_psi = np.pi/2
        

    def calc_drop_location(self,h,ptarget,v_des,into_wind = False):
        # Calculate Directional Velocity
        # V = [self.Va_hat * np.cos(self.psi_hat),
        #      self.Va_hat * np.sin(self.psi_hat),
        #      0]
        if into_wind == True:
            magV = np.linalg.norm(v_des)
            V = [magV * np.cos(self.wind_psi), -magV * np.sin(self.wind_psi), 0]
            self.approach_angle = -deepcopy(self.wind_psi)

        else:
            V = v_des # Desired Velocity to use in drop equations (i.e. how do you want to approach)
            self.approach_angle = np.arctan2(v_des[1],v_des[0])

        # Initial Position
        Z = [0,0,0]

        # Initial Time
        Tdrop = 0.0

        # Calculate Drop Position and Time
        h = -h
        while Z[2] < -h: 
            # Calculate dynamic forces on the object
            q = self.rho * self.S_object * self.Cd_object / (2 * self.mass_object)
            A = [-q * V[0]**2 * np.sign(V[0]),
                 -q * V[1]**2 * np.sign(V[1]),
                 self.g - q * V[2]**2 * np.sign(V[2])]

            # Solve Difference Equations for new velocity and position
            vprev = V
            V = np.add(V,np.multiply(A,self.dt))
            Z = np.add(Z,np.multiply(vprev,self.dt))

            # Update Time vector
            Tdrop += self.dt

            
        self.pdrop = np.add(np.add(ptarget, -Tdrop * self.wind_avg), np.multiply(-1,Z[0:2])) # Ground Frame?
        self.pdrop = [self.pdrop[0],self.pdrop[1],h]

        print 'Drop Location',self.pdrop

    def release_triggered(self,x,wind,target):
        # Determine True Location
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x

        R = [[np.cos(th) * np.cos(psi),np.sin(phi) * np.sin(th) * np.cos(psi) - np.cos(phi) * np.sin(psi),np.cos(phi) * np.sin(th) * np.cos(psi) + np.sin(phi) * np.sin(psi)],
             [np.cos(th) * np.sin(psi),np.sin(phi) * np.sin(th) * np.sin(psi) + np.cos(phi) * np.cos(psi),np.cos(phi) * np.sin(th) * np.sin(psi) - np.sin(phi) * np.cos(psi)],
             [-np.sin(th), np.sin(phi) * np.cos(th), np.cos(phi) * np.cos(th)]]
        vel_inertial = np.matmul(R,[u,v,w])
        u,v,w = vel_inertial
        V = [round(u,2),round(v,2),0.0]
        
        wn,we,wd,_,_,_ = wind
        Z = [0,0,0]
        Tdrop = 0.0
        q = self.rho * self.S_object * self.Cd_object / (2 * self.mass_object)

        p_hit = [pn,pe,-pd]
        while Z[2] < pd: 
            # Calculate dynamic forces on the object
            A = [-q * V[0]**2 * np.sign(V[0]),
                 -q * V[1]**2 * np.sign(V[1]),
                 self.g - q * V[2]**2 * np.sign(V[2])]


            # Solve Difference Equations for new velocity and position
            vprev = deepcopy(V)
            V = np.add(V,np.multiply(A,self.dt))
            Z = np.add(Z,np.multiply(vprev,self.dt))
            # state[0] += Z[0]
            # state[1] += Z[1]
            # state[2] += Z[2]
            # state[3] += V[0]
            # state[4] += V[1]
            # state[5] += V[2]

            # Update Time vector
            Tdrop += self.dt

            
        # Actual hit location
        p_hit += np.add(np.multiply(Tdrop,[wn,we,wd]),Z)
        self.hit = p_hit
        time.sleep(0.11)
        error = np.sqrt((target[0] - self.hit[0])**2 + (target[1] - self.hit[1])**2)
        print '\nAchieved:\t',self.hit,'\tError:\t',error,'\n'

    def bottle_eom(self,x,t):
        # Drag Forces
        q = self.rho * self.S_object * self.Cd_object / (2 * self.mass_object)

        # Equations of Motion
        xdot = [x[3],
                x[4],
                x[5],
                -q * x[3]**2,
                -q * x[4]**2,
                self.g -q * x[5]**2]

        return xdot

    def plan_path(self,p,t,chi):
        extend = [t[0] + 1000 * np.cos(self.approach_angle), t[1] + 1000 * np.sin(self.approach_angle),t[2],0]
        P = [[p[0],p[1],p[2],chi],[t[0],t[1],t[2],self.approach_angle],extend]
        return P

    def in_drop(self,p):
        d = np.sqrt((p[0] - self.pdrop[0])**2 + (p[1] - self.pdrop[1])**2)

        if d < 20.0:
            return True
        else:
            return False
        

    def plane_path_plot(self,w,p):
        pass
