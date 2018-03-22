from pathFollowing import *

class pathManager(pathFollow):
    def __init__(self,va0 = 30.0,model=1):
        pathFollow.__init__(self,va0,model)

        self.W = 0.0 # Waypoint path
        self.P = 0.0 # Dubins Configuration Path
        self.Rmin = 100.0 # Minimum turn radius

    def algorithm_5(self,W,p):
        '''
        Inputs: waypoints W
        MAV position p
        '''
        if len(W) < 3:
            print 'Not enough waypoints, setting r [0,0,0] and q [1,0,0]'
            r = [0,0,0]
            q_prev = [1,0,0]
        else:
            if W != self.W:
                self.W = W
                self.i = 1

            r = W[self.i - 1]
            q_prev = np.subtract(W[self.i],W[self.i - 1]) / np.linalg.norm(np.subtract(W[self.i],W[self.i - 1]))
            q = np.subtract(W[self.i + 1],W[self.i]) / np.linalg.norm(np.subtract(W[self.i + 1],W[self.i]))
            n = np.add(q,q_prev) / np.linalg.norm(np.add(q,q_prev))

            # if in half plane
            if self.half_plane(self.W[self.i],n,p) == True:
                self.i += 1

        return r, q_prev

    def algorithm_6(self,W,p,R):
        '''
        Inputs: waypoints W
        MAV position p
        fillet Radius R
        '''
        if len(W) < 3:
            print 'Not enough waypoints, setting r and q to default'
            r = [0,0,0]
            q = [1,0,0]
            flag = 1
            rho = R
            c = [0,0,0]
            lamb = 1
        else:
            if W != self.W:
                self.W = W
                self.i = 1
                self.path_state = 1

            q_prev = np.subtract(W[self.i],W[self.i - 1]) / np.linalg.norm(np.subtract(W[self.i],W[self.i - 1]))
            q = np.subtract(W[self.i + 1],W[self.i]) / np.linalg.norm(np.subtract(W[self.i + 1],W[self.i]))
            varrho = np.arccos(np.dot(-np.transpose(q_prev),q))

            if self.path_state == 1:
                flag = 1
                r = W[self.i - 1]
                q = q_prev
                z = np.subtract(W[self.i],np.multiply((R / np.tan(varrho/2)),q_prev))
                lamb = 1
                rho = R
                c = [0,0,0]

                if self.half_plane(z,q_prev,p) == True:
                    self.path_state = 2

            elif self.path_state == 2:
                flag = 2
                c = np.subtract(W[self.i],np.multiply((R / np.sin(varrho/2)),np.subtract(q_prev,q)/np.linalg.norm(np.subtract(q_prev,q))))
                rho = R
                lamb = np.sign(q_prev[0] * q[1] - q_prev[1] * q[0])
                z = np.add(W[self.i],np.multiply((R / np.tan(varrho/2)),q))
                r = [0,0,0]

                if self.half_plane(z,q,p) == True:
                    self.path_state = 1

        if flag == 1:
            return 'straight',r,q,c,rho,lamb
        elif flag == 2:
            return 'orbit',r,q,c,rho,lamb

    def algorithm_8(self,P,p,R):
        '''
        Inputs: configuration path P (waypoints with course angle)
        MAV position p
        Turn radius R
        '''
        if len(P) < 3:
            print 'Not enough waypoints, setting r and q to default'
            r = [0,0,0]
            q = [1,0,0]
            flag = 1
            rho = R
            c = [0,0,0]
            lamb = 1
        else:
            if P != self.P:
                self.W = W
                self.i = 1
                self.path_state = 1

        # Find dubins parameters
        L,cs,lambs,ce,lambe,z1,q1,z2,z3,q3 = self.find_dubins_parameters(P[self.i-1],P[self.i],R)

        if self.path_state == 1:
            flag = 2
            c = cs
            rho = R
            lamb = lambs

            if self.half_plane(z1,q1,p) == True:
                self.path_state = 2
                
        elif self.path_state == 2:
            if self.half_plane(z1,q1,p) == True:
                self.path_state = 3

        elif self.path_state == 3:
            flag = 1
            r = z1
            q = q1

            if self.half_plane(z2,q1,p) == True:
                self.path_state = 4

        elif self.path_state == 4:
            flag = 2
            c = ce
            rho = R
            lamb = lambe

            if self.half_plane(z3,-q3,p) == True:
                self.path_state = 5

        elif self.path_state == 5:
            if self.half_plane(z3,q3,p) == True:
                self.path_state = 1
                self.i += 1
                L,cs,lambs,ce,lambe,z1,q1,z2,z3,q3 = self.find_dubins_parameters(P[self.i-1],P[self.i],R)

        if flag == 1:
            return 'straight',r,q,c,rho,lamb
        elif flag == 2:
            return 'orbit',r,q,c,rho,lamb

    def find_dubins_parameters(prev,cur,R):
        ps = prev[0:3]
        pe = cur[0:3]
        chis = prev[3]
        chie = cur[3]
        if np.linalg.norm(np.subtract(ps,pe)) < 3 * R:
            return None
        elif R < self.Rmin:
            return None
        Rzpi = [[0,1,0],[-1,0,0],[0,0,1]]
        Rzpi_ = [[0,-1,0],[1,0,0],[0,0,1]]
        crs = np.add(ps,np.matmul(np.mutliply(R,Rzpi),[np.cos(chis),np.sin(chis),0]))
        cls = np.add(ps,np.matmul(np.mutliply(R,Rzpi_),[np.cos(chis),np.sin(chis),0]))
        cre = np.add(pe,np.matmul(np.mutliply(R,Rzpi),[np.cos(chie),np.sin(chie),0]))
        cle = np.add(pe,np.matmul(np.mutliply(R,Rzpi_),[np.cos(chie),np.sin(chie),0]))

        # Path 1
        vartheta = np.arccos(np.dot(crs,cre)/(np.linalg.norm(crs)*np.linalg.norm(cre)))
        L1 = np.linalg.norm(np.subtract(crs,cre)) + R * (2 * np.pi + vartheta - np.pi/2 - chis + np.pi/2) + R * (2 * np.pi + chie - np.pi/2 - vartheta + np.pi/2)

        # Path 2
        vartheta = np.arccos(np.dot(crs,cle)/(np.linalg.norm(crs)*np.linalg.norm(cle)))
        l = np.linalg.norm(np.subtract(cle,crs))
        vartheta2 = vartheta - np.pi/2 + np.arcsin(2 * R / l)
        L2 = np.sqrt(l**2 - 4 * R**2) + R * (2 * np.pi + vartheta2 - chis + np.pi/2) + R * (2 * np.pi + vartheta2 + np.pi - chie - np.pi/2)

        # Path 3
        vartheta = np.arccos(np.dot(cls,cre)/(np.linalg.norm(cls)*np.linalg.norm(cre)))
        l = np.linalg.norm(np.subtract(cre,cls))
        vartheta2 = np.arccos(2 * R / l)
        L3 = np.sqrt(l**2 - 4 * R**2) + R * (2 * np.pi + chis + np.pi/2 - vartheta - vartheta2) + R * (2 * np.pi + chie - np.pi/2 - vartheta - vartheta2 + np.pi)

        # Path 4
        vartheta = np.arccos(np.dot(cls,cle)/(np.linalg.norm(cls)*np.linalg.norm(cle)))
        L4 = np.linalg.norm(np.subtract(cls,cle)) + R * (2 * np.pi + chis + np.pi/2 - vartheta - np.pi/2) + R * (2 * pi + vartheta + np.pi/2 - chie - np.pi/2)

        L = np.min([L1,L2,L3,L4])
        arg = np.argmin([L1,L2,L3,L4])
        e1 = [1,0,0]

        if arg == 0:
            cs = crs
            lambs = 1
            ce = cre
            lambe = 1
            q1 = np.subtract(ce,cs)/np.linalg.norm(np.subtract(ce,cs))
            z1 = np.add(cs,np.matmul(np.multiply(R,Rzpi_),q1))
            z2 = np.add(ce,np.matmul(np.multiply(R,Rzpi_),q1))
        elif arg == 1:
            cs = crs
            lambs = 1
            ce = cle
            lambe = -1
            l = np.linalg.norm(np.subtract(ce,cs))
            vartheta = np.arccos(np.dot(ce,cs)/(np.linalg.norm(cs)*np.linalg.norm(ce)))
            vartheta2 = vartheta - np.pi/2 + np.arcsin(2*R/l)
            q1 = np.matmul(self.rotmatz(vartheta2 + np.pi/2),e1)
            z1 = np.add(cs,np.matmul(np.multiply(R,self.rotmatz(vartheta2)),e1))
            z2 = np.add(ce,np.matmul(np.multiply(R,self.rotmatz(vartheta2 + np.pi)),e1))
        elif arg == 2:
            cs = cls
            lambs = -1
            ce = cre
            lambe = +1
            l = np.linalg.norm(np.subtract(ce,cs))
            vartheta = np.arccos(np.dot(ce,cs)/(np.linalg.norm(cs)*np.linalg.norm(ce)))
            vartheta2 = vartheta - np.pi/2 + np.arcsin(2*R/l)
            q1 = np.matmul(self.rotmatz(vartheta + vartheta2 - np.pi/2),e1)
            z1 = np.add(cs,np.matmul(np.multiply(R,self.rotmatz(vartheta + vartheta2)),e1))
            z2 = np.add(ce,np.matmul(np.multiply(R,self.rotmatz(vartheta + vartheta2 - np.pi)),e1))
        elif arg == 3:
            cs = cls
            lambs = 1
            ce = cle
            lambe = 1
            q1 = np.subtract(ce,cs)/np.linalg.norm(np.subtract(ce,cs))
            z1 = np.add(cs,np.matmul(np.multiply(R,Rzpi),q1))
            z2 = np.add(ce,np.matmul(np.multiply(R,Rzpi),q1))

        z3 = pe
        q3 = np.matmul(self.rotmatz(chie),e1)

        return L,cs,lambs,ce,lambe,z1,q1,z2,z3,q3
            
    def rotmatz(self,angle):
        R = [[np.cos(angle),np.sin(angle),0],
             [-np.sin(angle),np.cos(angle),0],
             [0,0,1]]
        return R
            
        
    def half_plane(r,n,p):
        h = np.dot(np.transpose(np.subtract(p,r)),n)
        if h >= 0:
            return True
        else:
            return False
            
        
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

