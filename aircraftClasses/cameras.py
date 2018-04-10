import time
from pathPlanner import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class Cameras(pathPlanner):
    def __init__(self,va0 = 30.0,model=1):
        pathPlanner.__init__(self,va0,model)

        self.k_az = 1.0
        self.k_el = 1.0

    def point_gimbal(self,p_obj,p_mav,phi,th,psi,alpha_az,alpha_el):
        li_d = np.subtract(p_obj,p_mav)
        Rb_i = self.get_rotation_matrix([phi,th,psi])
        lb_d = 1/np.linalg.norm(li_d) * np.matmul(Rb_i,li_d)

        alphac_az = np.arctan(lb_d[1]/lb_d[0])
        alphac_el = np.arcsin(lb_d[2])

        u_az = self.k_az * (alphac_az - alpha_az)
        u_el = self.k_el * (alphac_el - alpha_el)

        return u_az,u_el,alphac_az,alphac_el

    def geolocation(self,p_obj,p_mav,phi,th,psi,az,el):
        Rcg = [[0,0,1],
               [1,0,0],
               [0,1,0]]
        Rbi = self.get_rotation_matrix([phi,th,psi])
        u_az,u_el,azc,elc = self.point_gimbal(p_obj,p_mav,phi,th,psi,az,el)
        Rgb = self.rot_b_to_g(azc,elc)

        varrho = np.dot([0,0,1],li)
        L = -p_mav[2] / varrho

    def rot_b_to_g(self,az,el):
        Rbg1 = [[np.cos(az),np.sin(az),0],
                [-np.sin(az),np.cos(az),0],
                [0,0,1]]
        Rg1g = [[np.cos(el),0,-np.sin(el)],
                [0,1,0],
                [np.sin(el),0,np.cos(el)]]
        return np.matmul(Rg1g,Rbg1)
        

        

        
