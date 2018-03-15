from kalmanFilter import *
import pyqtgraph as pg

# class guidanceModel(kalmanFilter):
#     def __init__(self,x0,trim):
#         kalmanFilter.__init__(self,x0,trim)
class guidanceModel():
    def __init__(self,va0 = 30.0,model=1):
        # Initial States
        if model == 1:
            self.x = [0,0,0,0,0,0,va0] # pn,pe,chi,chidot,h,hdot,va
        elif model == 2:
            self.x = [0,0,0,0,0,va0,0] # pn,pe,psi,h,hdot,va,phi
        self.dt = 0.01
        self.t_sim = 0.0
        self.va0 = va0
        self.first = True

        # Parameters to tune
        self.bchidot = 4.5
        self.bchi = 5.0
        self.bhdot = 0.8
        self.bh = 1.9
        self.bva = 6.5
        self.bphi = 3.0
        
    def model1(self,x,t,wn,we,wd,chidot_c,chi_c,hdot_c,h_c,va_c):
        # Rename Inputs
        pn,pe,chi,chidot,h,hdot,va = x

        # Calculate psi (Assuming gamma a = 0)
        psi = chi - np.arcsin(1/va * (-wn * np.sin(chi) + we * np.cos(chi)))

        # Eq 9.19
        pndot = va * np.cos(psi) + wn
        pedot = va * np.sin(psi) + we
        chiddot = self.bchidot * (chidot_c - chidot) + self.bchi * (chi_c - chi)
        hddot = self.bhdot * (hdot_c - hdot) + self.bh * (h_c - h)
        vadot = self.bva * (va_c - va)

        return [pndot,pedot,chidot,chiddot,hdot,hddot,vadot]

    def model2(self,x,t,wn,we,wd,hdot_c,h_c,va_c,phi_c):
        # Rename Inputs
        pn, pe, psi, h, hdot, va, phi = x

        # Eq 9.20
        pndot = va * np.cos(psi) + wn
        pedot = va * np.sin(psi) + we
        psidot = 9.81 / va * np.tan(phi)
        hddot = self.bhdot * (hdot_c - hdot) + self.bh * (h_c - h)
        vadot = self.bva * (va_c - va)
        phidot = self.bphi * (phi_c - phi)

        # Turn Radius
        R = va**2 / (9.81 * np.tan(phi))
        print R

        return [pndot,pedot,psidot,hdot,hddot,vadot,phidot]

    def simulate_states(self,winds,commands,model=1):
        # Wind Inputs
        wn,we,wd = winds

        # Check which model and solve ode
        if model == 1:
            chidot_c, chi_c, hdot_c, h_c, va_c, _ = commands
            sol = odeint(self.model1,self.x,[0.0,self.dt],
                         args=(wn,we,wd,chidot_c,chi_c,hdot_c,h_c,va_c))
        elif model == 2:
            _, _, hdot_c, h_c, va_c, phi_c = commands
            sol = odeint(self.model2,self.x,[0.0,self.dt],
                         args=(wn,we,wd,hdot_c,h_c,va_c,phi_c))

        # Update values
        # Model 1 states: pn,pe,chi,chidot,h,hdot,va
        # Model 2 states: pn,pe,psi,h,hdot,va,phi
        self.x = sol[-1]


    def plot(self,cmd,state,auto,labels):
        if self.first == True:
            # Setup Window
            self.app = pg.QtGui.QApplication([])
            self.app.aboutToQuit.connect(self.stop)
            self.win = pg.GraphicsWindow(size=(1200,400))
            self.win.setWindowTitle('States')
            self.win.setInteractive(True)

            # Plots and Subplots
            self.p = []
            for i in xrange(3):
                self.p.append(self.win.addPlot())
                self.p[i].setLabel('bottom','Time',units='s')
                self.p[i].setLabel('left',labels[i])
                self.p[i].addLegend()

            # Curves
            self.commandcurve = []
            self.statecurve = []
            self.autocurve = []

            # Data
            self.command = [[],[],[]]
            self.state = [[],[],[]]
            self.autostate = [[],[],[]]

            # Pud data on curves
            for i in xrange(3):
                pen = pg.mkPen(color=(0,255,0),style=pg.QtCore.Qt.DashLine)
                self.commandcurve.append(self.p[i].plot(pen=pen,name='Command'))
                pen = pg.mkPen(color=(0,0,255))
                self.statecurve.append(self.p[i].plot(pen=pen,name='Ch. 9 State'))
                pen = pg.mkPen(color=(255,0,0))
                self.autocurve.append(self.p[i].plot(pen=pen,name='AutoPilot'))
                self.command[i].append(cmd[i])
                self.state[i].append(state[i])
                self.autostate[i].append(auto[i])

            self.first = False
        else:
            time = np.linspace(0,(len(self.command[0]) + 1) * self.dt, len(self.command[0]) + 1)
            for i in xrange(len(cmd)):
                self.command[i].append(cmd[i])
                self.state[i].append(state[i])
                self.autostate[i].append(auto[i])
                self.commandcurve[i].setData(time,self.command[i])
                self.statecurve[i].setData(time,self.state[i])
                self.autocurve[i].setData(time,self.autostate[i])

            # Update Plot
            self.app.processEvents()

            
    def stop(self):
        sys.exit()
