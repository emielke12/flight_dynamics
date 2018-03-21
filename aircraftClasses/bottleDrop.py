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

    def wind_estimate(self,x0):
        # Estimate wind from current state or state estimate (if using kalman filter)
        wn = self.Vg_hat * np.cos(self.chi_hat) - self.Va_hat * np.cos(self.psi_hat)
        we = self.Vg_hat * np.sin(self.chi_hat) - self.Va_hat * np.sin(self.psi_hat)

        if len(self.wind_array) < 10:
            self.wind_array.append([wn,we])
        else:
            self.wind_array = np.roll(self.wind_array,-1,axis=0)
            self.wind_array[-1] = [wn,we]
        self.wind_avg = np.average(self.wind_array,axis=0)
        self.wind_psi = np.arctan2(self.wind_avg[1],self.wind_avg[0])

    def calc_drop_location(self,h,ptarget,v_des):
        # Calculate Directional Velocity
        # V = [self.Va_hat * np.cos(self.psi_hat),
        #      self.Va_hat * np.sin(self.psi_hat),
        #      0]
        V = v_des # Desired Velocity to use in drop equations (i.e. how do you want to approach)

        # Initial Position
        Z = [0,0,0]

        # Initial Time
        Tdrop = 0.0

        # Calculate Drop Position and Time
        while Z[2] < h: 
            # Calculate dynamic forces on the object
            q = self.rho * self.S_object * self.Cd_object / (2 * self.mass_object)
            A = [-q * V[0]**2,
                 -q * V[1]**2,
                 self.g - q * V[2]**2]

            # Solve Difference Equations for new velocity and position
            vprev = V
            V = np.add(V,np.multiply(A,self.dt))
            Z = np.add(Z,np.multiply(vprev,self.dt))

            # Update Time vector
            Tdrop += self.dt
            
#         # Wind direction
#         Rpsi = [[np.cos(self.wind_psi), -np.sin(self.wind_psi), 0],
#                 [np.sin(self.wind_psi), np.cos(self.wind_psi), 0],
#                 [0, 0, 1]]

        # Update drop location
#         pdrop = np.add(np.add(ptarget, -Tdrop * self.wind_avg), -np.matmul(Rpsi,Z)[0:2])
        self.pdrop = np.add(np.add(ptarget, -Tdrop * self.wind_avg), -Z[0:2]) # Ground Frame?
        # print 'Drop Location',self.pdrop

    def release_triggered(self,x,wind):
        # Determine True Location
        pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
        wn,we,wd,_,_,_ = wind
        self.true_bottle_path = [[pn],[pe],[-pd]]
        state = [pn,
                 pe,
                 -pd,
                 u - wn,
                 v - we,
                 w - wd]
        time = 0.0
        while state[2] > 0: 
            sol = odeint(self.bottle_eom,state,[0.0, self.dt])
            time += self.dt
            state = [sol[-1][0] + self.dt * wn,
                     sol[-1][1] + self.dt * we,
                     sol[-1][2] + self.dt * wd,
                     sol[-1][3],
                     sol[-1][4],
                     sol[-1][5]]
            for i in xrange(3):
                self.true_bottle_path[i].append(state[i])
    
        # Actual hit location
        self.hit = [self.true_bottle_path[0][-1],self.true_bottle_path[1][-1]]
        print '\nAchieved Hit at: ',self.hit,'\n'

    def bottle_eom(self,x,t):
        # Drag Forces
        q = self.rho * self.S_object * self.Cd_object / (2 * self.mass_object)

        # Equations of Motion
        xdot = [x[3],
                x[4],
                -x[5],
                -q * x[3]**2,
                -q * x[4]**2,
                self.g -q * x[5]**2]

        return xdot

    def in_drop(self):
        d = np.sqrt((self.x_hat[0] - self.pdrop[0])**2 + (self.x_hat[1] - self.pdrop[1])**2)
        if d < 5.0:
            return True
        else:
            return False
        
