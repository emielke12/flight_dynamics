from guidance import *

class pathFollow(guidanceModel):
    def __init__(self,va0 = 30.0,model=1):
        guidanceModel.__init__(self,va0,model)

        # Gains
        self.chi_inf =  60 * np.pi / 180.0
        self.k_path = 0.05
        self.k_orbit = 5.0
        self.lamb = 1 # 1 is clockwise, -1 is counter clockwise

    def follow_straight_line(self,r,q,p,chi):
        ''' Inputs:
        r -- flight origin
        q -- unit vector of path direction
        p -- MAV position
        chi -- course angle
        '''
        # Path Course Angle
        chi_q = np.arctan2(q[1],q[0])
        while chi_q - chi < -np.pi:
            chi_q += 2 * np.pi
        while chi_q - chi > np.pi:
            chi_q -= 2 * np.pi

        # Transformation from inertial frame to path frame
        R_p = [[np.cos(chi_q), np.sin(chi_q), 0],
               [-np.sin(chi_q), np.cos(chi_q), 0],
               [0, 0, 1]]

        # Error in inertial frame
        e_pi = [p[0] - r[0],
                p[1] - r[1],
                p[2] - r[2]]

        # k vector
        k = [0,0,1]

        # Vector normal to q-k plane
        n = np.cross(q,k) / np.linalg.norm(np.cross(q,k))

        # Projection of e_p onto q-k plane
        s_i = np.subtract(e_pi, np.multiply(np.dot(e_pi,n),n))

        # Desired altitude
        h_c = -r[2] - np.sqrt(s_i[0]**2 + s_i[1]**2) * (q[2] / (np.sqrt(q[0]**2 + q[1]**2)))

        # Error in path frame (y-direction)
        e_p = np.matmul(R_p,e_pi)
        e_py = e_p[1]

        # Desired course
        chi_d = -self.chi_inf * 2 / np.pi * np.arctan(self.k_path * e_py)

        # Commanded course angle
        chi_c = chi_q + chi_d

        return h_c, chi_c

    def follow_orbit(self,c,rho,p,chi):
        ''' Inputs:
        c -- orbit center
        rho -- orbit radius
        lamb -- orbit direction
        p -- MAV position
        chi -- course angle
        '''
        # Desired altitude
        h_c = c[2]

        # Radial distance from center of desired orbit
        d = np.sqrt((p[0] - c[0])**2 + (p[1] - c[1])**2)

        # phase angle of relative position (from north)
        varphi = np.arctan2(p[1] - c[1], p[0] - c[0])

        # Checks
        while varphi - chi < -np.pi:
            varphi += 2 * np.pi
        while varphi - chi > np.pi:
            varphi -= 2 * np.pi    

        # Commanded course angle
        chi_c =  varphi + self.lamb * (np.pi / 2 + np.arctan(self.k_orbit * (d - rho) / rho))
        
        return h_c, chi_c
        
    def path_plot(self,r,q,c,rho,p,path_type='straight'):
        if self.t_sim < 2 * self.dt:
            # Setup Window
            self.path_win = pg.GraphicsWindow(size=(800,800))
            self.path_win.setWindowTitle('Path')
            self.path_win.setInteractive(True)

            # Plots
            self.path_p = self.path_win.addPlot()
            self.path_p.addLegend()

            # Data
            self.achieved_path = [[],[]]

            # Put data on curves
            pen = pg.mkPen(color=(0,255,0),style=pg.QtCore.Qt.DashLine)
            self.command_path_curve = self.path_p.plot(pen=pen,name='Command')
            pen = pg.mkPen(color=(0,0,255))
            self.achieved_path_curve = self.path_p.plot(pen=pen,name='Achieved')
            
            self.achieved_path[0].append(p[0])
            self.achieved_path[1].append(p[1])
            chi_q = np.arctan2(q[1],q[0])
            
        else:
            self.achieved_path[0].append(p[0])
            self.achieved_path[1].append(p[1])
            self.achieved_path_curve.setData(self.achieved_path[1],self.achieved_path[0])
            
            chi_q = np.arctan2(q[1],q[0])

            if path_type == 'straight':
                self.command_path = [[r[0],r[0] + self.t_sim * self.x[-1] * np.cos(chi_q)],
                                     [r[1],r[1] + self.t_sim * self.x[-1] * np.sin(chi_q)]]
                self.command_path_curve.setData(self.command_path[1],self.command_path[0])
            elif path_type == 'orbit':
                theta = np.linspace(0,2 * np.pi, 50)
                self.command_path = [np.add(c[0],np.multiply(rho,np.cos(theta))),
                                     np.add(c[1],np.multiply(rho,np.sin(theta)))]
                self.command_path_curve.setData(self.command_path[1],self.command_path[0])
            
            # Update Plot
            self.app.processEvents()

