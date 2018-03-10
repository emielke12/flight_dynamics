from kalmanFilter import *
import pyqtgraph as pg

class bottleDrop(kalmanFilter):
    def __init__(self,x0,trim):
        kalmanFilter.__init__(self,x0,trim)
        self.wind_array = []
        self.wind_avg = [0.0,0.0]

        # If object is sphere and low speed
        self.radius = 0.25 # m
        self.mass_object = 0.5 # kg
        self.S_object = np.pi * self.radius ** 2
        self.Cd_object = 0.5

    def wind_estimate(self,x0):
        # Estimate wind from current state or state estimate (if using kalman filter)
        wn = self.Vg_hat * np.cos(self.chi_hat) - self.Va_hat * np.cos(self.psi_hat)
        we = self.Vg_hat * np.sin(self.chi_hat) - self.Va_hat * np.sin(self.psi_hat)

        if len(self.wind) < 10:
            self.wind_array.append([wn,we])
        else:
            self.wind_array = np.roll(self.wind_array,-1,axis=0)
            self.wind_array[-1] = [wn,we]
        self.wind_avg = np.average(self.wind_array,axis=0)

    def calc_drop_location(self,h,ptarget,x,wind):
        # Calculate Directional Velocity
        V = [self.Va_hat * np.cos(self.psi_hat),
             self.Va_hat * np.sin(self.psi_hat),
             0]

        # Initial Position
        Z = [0,0,0]

        # Initial Time
        Tdrop = 0.0

        # Drop Location
        pdrop = [0,0]

        while Z[2] < h: 
            # Calculate dynamic forces on the object
            q = self.rho * self.S_object * self.Cd_object / (2 * self.mass_object)
            A = [-q * V[0]**2,
                 -q * V[1]**2,
                  self.g - q * V[2]**2]

            # Solve Difference Equations for new velocity and position
            V = np.add(V,A * self.dt)
            Z = np.add(Z,V * self.dt)

            # Update Time vector
            Tdrop += self.dt
            
        # Wind direction
        psi_wind = np.atan2(self.wind_avg[1],self.wind_avg[0])
        Rpsi = [[np.cos(psi_wind), -np.sin(psi_wind), 0],
                [np.sin(psi_wind), np.cos(psi_wind), 0],
                [0, 0, 1]]

        # Update drop location
        pdrop = np.add(np.add(ptarget, -Tdrop * self.wind_avg), -np.matmul(Rpsi,Z)[0:2])
        print 'Drop Location',pdrop

        # Something is wrong with this stuff
#         # Determine True Location
#         self.true_bottle_path = [pdrop[0],pdrop[1],-x[2]]
#         while self.true_bottle_path[2] < h: 
#             sol = odeint(self.bottle_eom,x,[0.0, self.dt],args=(wind[0],wind[1],wind[2]))
#             for i in xrange(3):
#                 x[i + 3] = sol[-1][i + 3]
#                 self.true_bottle_path[i].append(sol[-1][i])
    
#         # Actual hit location
#         self.hit = [self.true_bottle_path[0][-1],self.true_bottle_path[1][-1]]
#         print '\nAchieved Hit at: ',self.hit,'\n'

#     def bottle_eom(self,x,wn,we,wd):
#         # Rename states
#         pn,pe,pd,u,v,w,phi,th,psi,p,q,r = x
#         wn,we,wd = wind

#         # Uses True States
#         v = [u - wn,
#              v - we,
#              w - wd]
        
#         # Drag Forces
#         q = self.rho * self.S_object * self.Cd_object / (2 * self.mass_object)

#         xdot = [v[0],
#                 v[1],
#                 v[2],
#                 -q * v[0]**2
#                 -q * v[1]**2
#                 self.g -q * v[2]**2]

#         return xdot
