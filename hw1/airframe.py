from mpl_toolkits.mplot3d import Axes3D as a3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection as pc
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

class AirCraftDrawing():
    def __init__(self):
        # Aircraft Parameters
        self.fuse_l1 = 1.0
        self.fuse_l2 = 0.5
        self.fuse_l3 = 3.0
        self.wing_l = 0.5
        self.fuse_h = 0.5
        self.fuse_w = 0.5
        self.wing_l = 1.0
        self.wing_w = 3.0
        self.tailwing_l = 0.5
        self.tailwing_w = 1.5
        self.tail_h = 0.5

        # Figure
        plt.close('all')
        self.fig = plt.figure()
        self.ax = a3(self.fig)

    def create_wing(self,R,t):
        x = [0, 0, -self.wing_l, -self.wing_l]
        y = [self.wing_w/2, -self.wing_w/2, -self.wing_w/2, self.wing_w/2]
        z = [0,0,0,0]
        new = np.split(np.dot(R,np.column_stack((np.column_stack((x,y)),z)).T),3)
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1]) 
        z = np.add(new[2][0],t[2])
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts))

    def create_fuselage(self,R,t):
        x = [-self.fuse_l3, self.fuse_l2, self.fuse_l2,-self.fuse_l3, self.fuse_l2, self.fuse_l2,-self.fuse_l3, self.fuse_l2, self.fuse_l2,-self.fuse_l3, self.fuse_l2, self.fuse_l2]
        y = [0, self.fuse_w/2, -self.fuse_w/2, 0, self.fuse_w/2, -self.fuse_w/2, 0, self.fuse_w/2, self.fuse_w/2, 0, -self.fuse_w/2, -self.fuse_w/2]
        z = [0, self.fuse_h/2, self.fuse_h/2, 0, -self.fuse_h/2, -self.fuse_h/2, 0, -self.fuse_h/2, self.fuse_h/2, 0, -self.fuse_h/2, self.fuse_h/2]
        new = np.split(np.dot(R,np.column_stack((np.column_stack((x,y)),z)).T),3)
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1]) 
        z = np.add(new[2][0],t[2])
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts))

    def create_nose(self,R,t):
        x = [self.fuse_l2,self.fuse_l1,self.fuse_l2,self.fuse_l2,self.fuse_l1,self.fuse_l2,self.fuse_l2,self.fuse_l1,self.fuse_l2,self.fuse_l2,self.fuse_l1,self.fuse_l2]
        y = [self.fuse_w/2,0,-self.fuse_w/2,self.fuse_w/2,0,self.fuse_w/2,-self.fuse_w/2,0,-self.fuse_w/2,self.fuse_w/2,0,-self.fuse_w/2]
        z = [self.fuse_h/2,0,self.fuse_h/2,self.fuse_h/2,0,-self.fuse_h/2,self.fuse_h/2,0,-self.fuse_h/2,-self.fuse_h/2,0,-self.fuse_h/2]
        new = np.split(np.dot(R,np.column_stack((np.column_stack((x,y)),z)).T),3)
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1]) 
        z = np.add(new[2][0],t[2])
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts))

    def create_tailwing(self,R,t):
        x = [-self.fuse_l3, -self.fuse_l3, -self.fuse_l3 + self.tailwing_l, -self.fuse_l3 + self.tailwing_l]
        y = [self.tailwing_w/2, -self.tailwing_w/2, -self.tailwing_w/2, self.tailwing_w/2]
        z = [0,0,0,0]
        new = np.split(np.dot(R,np.column_stack((np.column_stack((x,y)),z)).T),3)
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1]) 
        z = np.add(new[2][0],t[2])
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts))

    def create_tailfin(self,R,t):
        x = [-self.fuse_l3,-self.fuse_l3,-self.fuse_l3 + self.tailwing_l]
        y = [0,0,0]
        z = [0,self.tail_h,0]
        new = np.split(np.dot(R,np.column_stack((np.column_stack((x,y)),z)).T),3)
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1]) 
        z = np.add(new[2][0],t[2])
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts))


    def create_plane(self,eul = [0,0,0],t = [0,0,0]):
        R = self.get_rotation_matrix(eul)
        self.ax.clear()
        self.create_wing(R,t)
        self.create_fuselage(R,t)
        self.create_nose(R,t)
        self.create_tailwing(R,t)
        self.create_tailfin(R,t)

        plt.draw()
        
    def get_rotation_matrix(self,eul):
        psi = eul[0]
        th = eul[1]
        phi = eul[2]
        R = [[np.cos(th)*np.cos(psi), np.cos(th)*np.sin(psi), -np.sin(th)],
             [np.sin(phi)*np.sin(th)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(th)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(th)],
             [np.cos(phi)*np.sin(th)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(th)]]
        return R

    def animate(self,i):
        eul = [0,0,ang[i]]
        t = [0,0,0]
        self.create_plane(eul,t)
        print i


A = AirCraftDrawing()
ang = np.linspace(0,np.pi/2,100)
for j in xrange(90):
    ani = animation.FuncAnimation(A.fig,A.animate,len(ang),interval=10,repeat = False)
