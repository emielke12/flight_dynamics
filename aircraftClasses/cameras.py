import time
from pathPlanner import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class Cameras(pathPlanner):
    def __init__(self,va0 = 30.0,model=1):
        pathPlanner.__init__(self,va0,model)

        self.k_az = 10
        self.k_el = 10
        self.az0 = 0
        self.el0 = -np.pi/2
        self.cam_pix = 480
        self.cam_fov = 10 * np.pi / 180.0
        self.cam_f = 0.5 * self.cam_pix / np.tan(self.cam_fov/2)
        self.target_size = 2.0
        self.az_lim = 180 * np.pi / 180.0
        self.el_lim = 180 * np.pi / 180.0

    def point_gimbal(self,p_obj,p_mav,phi,th,psi,lcd,Rcg,Rgb):
        li_d = np.subtract(p_obj,p_mav)
        Rb_i = self.get_rotation_matrix(phi,th,psi)
        lb_d = 1/np.linalg.norm(li_d) * np.matmul(Rb_i,li_d)

        # RR = np.matmul(Rgb,Rcg)
        # lb_d = np.matmul(RR,lcd)

        # az_d = np.arctan2(lb_d[1]/lb_d[0],1)
        az_d = np.arctan(lb_d[1]/lb_d[0])
        el_d = np.arcsin(lb_d[2])
        
        u_az = self.k_az * (az_d - deepcopy(self.az0))
        u_el = self.k_el * (el_d - deepcopy(self.el0))

        return u_az,u_el

    def geolocation(self,p_obj,p_mav,phi,th,psi):
        Rcg = [[0,0,1],
               [1,0,0],
               [0,1,0]]
        Rbi = self.get_rotation_matrix(phi,th,psi)
        Rgb = self.rot_b_to_g(deepcopy(self.az0),deepcopy(self.el0))
        RR = np.matmul(np.matmul(Rbi,Rgb),Rcg)
        epsx,epsy,epss = self.camProj(p_obj,p_mav,Rbi,Rgb)
        lc_d = np.multiply(1 / np.sqrt(self.cam_f**2 + epsx**2 + epsy**2),[epsx,epsy,self.cam_f])
        u_az,u_el = self.point_gimbal(p_obj,p_mav,phi,th,psi,lc_d,Rcg,Rgb)

        self.determineNewAlpha(u_az,u_el)
        self.drawCamera(epsx,epsy)

        calc_pos = self.determine_location(p_mav,Rcg,Rbi,Rgb,lc_d)
        self.geo_error = calc_pos - p_obj

    def get_rotation_matrix(self,phi,th,psi):
        R = [[np.cos(th) * np.cos(psi),np.sin(phi) * np.sin(th) * np.cos(psi) - np.cos(phi) * np.sin(psi),np.cos(phi) * np.sin(th) * np.cos(psi) + np.sin(phi) * np.sin(psi)],
             [np.cos(th) * np.sin(psi),np.sin(phi) * np.sin(th) * np.sin(psi) + np.cos(phi) * np.cos(psi),np.cos(phi) * np.sin(th) * np.sin(psi) - np.sin(phi) * np.cos(psi)],
             [-np.sin(th), np.sin(phi) * np.cos(th), np.cos(phi) * np.cos(th)]]

        # Rvv1 = [[np.cos(psi), np.sin(psi), 0],
        #         [-np.sin(psi), np.cos(psi),0],
        #         [0,0,1]]
        # Rv1v2 = [[np.cos(th), 0, -np.sin(th)],
        #          [0,1,0],
        #          [np.sin(th), 0, np.cos(th)]]
        # Rv2b = [[1, 0, 0],
        #         [0, np.cos(phi), np.sin(phi)],
        #         [0, -np.sin(phi), np.cos(phi)]]
        # R = np.matmul(np.matmul(Rv2b,Rv1v2),Rvv1)
        return R


    def rot_b_to_g(self,az,el):
        Rbg1 = [[np.cos(az),np.sin(az),0],
                [-np.sin(az),np.cos(az),0],
                [0,0,1]]
        Rg1g = [[np.cos(el),0,-np.sin(el)],
                [0,1,0],
                [np.sin(el),0,np.cos(el)]]
        return np.matmul(Rg1g,Rbg1)
        

    def camProj(self,p_obj,p_mav,Rbi,Rgb):
        Rgc = [[0,1,0],
               [0,0,1],
               [1,0,0]]
        RR = np.matmul(np.matmul(Rgc,Rgb),Rbi)
        p_obj_v = np.subtract(p_obj,p_mav)
        p_obj_c = np.matmul(RR,p_obj_v)


        epsx = self.cam_f * p_obj_c[0] / p_obj_c[2]
        epsy = self.cam_f * p_obj_c[1] / p_obj_c[2]
        epss = self.cam_f * self.target_size / p_obj_c[2]

        tmp = self.cam_pix / 2 + epss
        if epsx < -tmp or epsx > tmp or epsy < -tmp or epsy > tmp:
            epsx = -9999
            epsy = -9999
            epss = 0

        return epsx,epsy,epss
                             
    def drawCamera(self,epsx,epsy):
        if self.t_sim < 2 * self.dt:
            self.cam_win = pg.GraphicsWindow(size=(800,800))
            self.cam_win.setWindowTitle('Camera')
            self.cam_win.setInteractive(True)

            self.cam_plot = self.cam_win.addPlot()
            self.cam_plot.setXRange(-self.cam_pix/2,self.cam_pix/2,padding=0)
            self.cam_plot.setYRange(-self.cam_pix/2,self.cam_pix/2,padding=0)
            self.obj_arrow = pg.ArrowItem(angle=0)
            self.cam_plot.addItem(self.obj_arrow)
            self.obj_arrow.setPos(epsy,epsx)

        else:
            self.obj_arrow.setPos(epsy,epsx)
            self.obj_arrow.setRotation(90)
            self.app.processEvents()
            
    def determineNewAlpha(self,uaz,uel):
        self.az0 = deepcopy(self.az0) + uaz * self.dt
        self.el0 = deepcopy(self.el0) + uel * self.dt

        if self.az0 < -self.az_lim:
            self.az0 = -self.az_lim
        elif self.az0 > self.az_lim:
            self.az0 = self.az_lim
            
        if self.el0 < -self.el_lim:
            self.el0 = -self.el_lim
        elif self.el0 > self.el_lim:
            self.el0 = self.el_lim

    def determine_location(self,p_mav,Rcg,Rbi,Rgb,lc):
        RR = np.matmul(np.matmul(np.transpose(Rbi),np.transpose(Rgb)),Rcg)
        righttop = np.matmul(RR,lc)
        rightbottom = np.dot([0,0,1],np.matmul(RR,lc))
        p_obj = np.add(p_mav, np.multiply(-p_mav[2],np.divide(righttop,rightbottom)))
        return p_obj

        
