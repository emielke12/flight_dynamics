from mpl_toolkits.mplot3d import Axes3D as a3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection as pc
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from Tkinter import *
import matplotlib.animation as animation
import numpy as np
import time,sys
from copy import deepcopy

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
        self.axes_length = 3.5

        # Make Figure
        self.bounds = 30
        self.bound_list = [10.0,10.0,10.0]
        self.make_figure()
        
        # Animation Values
        self.n_frames = 400
        self.interv = 10
        
        # # Visualization start points
        # self.phi = 0.0
        # self.th = 0.0
        # self.psi = 0.0
        # self.x = 0.0
        # self.y = 0.0
        # self.z = 0.0
        # self.init_phi = 0.0
        # self.init_th = 0.0
        # self.init_psi = 0.0
        # self.init_x = 0.0
        # self.init_y = 0.0
        # self.init_z = 0.0
        
    def make_figure(self,real_time = True):
        # Figure
        self.fig = plt.figure(figsize=(10,10)) # Figure Window size
        self.ax = a3(self.fig) # 3d axis
        self.ax.set_xbound(-self.bounds,self.bounds) # Axis bounds
        self.ax.set_ybound(-self.bounds,self.bounds)
        self.ax.set_zbound(-self.bounds,self.bounds)
        self.ax.view_init(elev=30,azim=140) # Axis view
        self.ax.dist=10 
        self.ax.invert_zaxis() # North east down view
        self.ax.invert_yaxis() # North east down view
        if real_time == True:
            self.fig.show()
            # self.root = Tk()
            # ele = Scale(self.root,from_=-40.0,to=40.0,orient=HORIZONTAL,label='Elevator',resolution=0.01)
            # ele.pack()
            # ail = Scale(self.root,from_=-40.0, to=40.0,orient=HORIZONTAL,label='Aileron',resolution=0.01)
            # ail.pack()
            # rud = Scale(self.root,from_=-40.0, to=40.0,orient=HORIZONTAL,label='Rudder',resolution=0.01)
            # rud.pack()
            # thr = Scale(self.root,from_=0.0, to=3.0,orient=HORIZONTAL,label='Thrust',resolution=0.01)
            # thr.pack()
            # self.sliders = [ele,ail,rud,thr]
                                                    


    def draw_line(self,px,py,pz,R,t,colour='r',lw = '5'):
        ''' Draw a line for axes'''
        # Rotate
        new = np.split(np.dot(R,np.column_stack((np.column_stack((px,py)),pz)).T),3)

        # Translate
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1])
        z = np.add(new[2][0],t[2])

        # Add to collection
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts,edgecolors=colour,linewidth=lw))
        
    def draw_2dpolygon(self,px,py,pz,R,t,colour='#929591'):
        '''Draw a polygon with inputs px,py,pz which are a list of vertices
        Rotate the points according to the specified rotation matrix R
        Translate the points according to the translation vector t
        Give polygon a color according to 'r' '''
        # Rotate
        new = np.split(np.dot(R,np.column_stack((np.column_stack((px,py)),pz)).T),3)

        # Translate
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1])
        z = np.add(new[2][0],t[2])

        # Add to collection
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts,facecolors=colour,edgecolors='k'))

    # These all create the different parts of the aircraft
    def create_wing(self,R,t):
        x = [0, 0, -self.wing_l, -self.wing_l]
        y = [-self.wing_w/2, self.wing_w/2, self.wing_w/2, -self.wing_w/2]
        z = [0,0,0,0]
        self.draw_2dpolygon(x,y,z,R,t)

    def create_fuselage(self,R,t):
        x = [[-self.fuse_l3, self.fuse_l2, self.fuse_l2],[-self.fuse_l3, self.fuse_l2, self.fuse_l2],[-self.fuse_l3, self.fuse_l2, self.fuse_l2],[-self.fuse_l3, self.fuse_l2, self.fuse_l2]]
        y = [[0, -self.fuse_w/2, self.fuse_w/2], [0, self.fuse_w/2, self.fuse_w/2], [0, self.fuse_w/2, -self.fuse_w/2], [0, -self.fuse_w/2, -self.fuse_w/2]]
        z = [[0, -self.fuse_h/2, -self.fuse_h/2], [0,-self.fuse_h/2, self.fuse_h/2], [0,self.fuse_h/2, self.fuse_h/2], [0,self.fuse_h/2, -self.fuse_h/2]]
        for i in xrange(len(x)):
            self.draw_2dpolygon(x[i],y[i],z[i],R,t)

    def create_nose(self,R,t):
        x = [[self.fuse_l2,self.fuse_l1,self.fuse_l2],[self.fuse_l2,self.fuse_l1,self.fuse_l2],[self.fuse_l2,self.fuse_l1,self.fuse_l2],[self.fuse_l2,self.fuse_l1,self.fuse_l2]]
        y = [[-self.fuse_w/2,0,self.fuse_w/2],[-self.fuse_w/2,0,-self.fuse_w/2],[self.fuse_w/2,0,self.fuse_w/2],[-self.fuse_w/2,0,self.fuse_w/2]]
        z = [[-self.fuse_h/2,0,-self.fuse_h/2],[-self.fuse_h/2,0,self.fuse_h/2],[-self.fuse_h/2,0,self.fuse_h/2],[self.fuse_h/2,0,self.fuse_h/2]]
        for i in xrange(len(x)):
            self.draw_2dpolygon(x[i],y[i],z[i],R,t)
        
    def create_tailwing(self,R,t):
        x = [-self.fuse_l3, -self.fuse_l3, -self.fuse_l3 + self.tailwing_l, -self.fuse_l3 + self.tailwing_l]
        y = [-self.tailwing_w/2, self.tailwing_w/2, self.tailwing_w/2, -self.tailwing_w/2]
        z = [0,0,0,0]
        self.draw_2dpolygon(x,y,z,R,t)
        
    def create_tailfin(self,R,t):
        x = [-self.fuse_l3,-self.fuse_l3,-self.fuse_l3 + self.tailwing_l]
        y = [0,0,0]
        z = [0,-self.tail_h,0]
        self.draw_2dpolygon(x,y,z,R,t)

    def create_axes(self,R,t):
        x = [0,self.axes_length]
        y = [0,0]
        z = [0,0]

        # Body Frame
        self.draw_line(x,y,z,R,t,'b',lw='3')
        self.draw_line(y,x,z,R,t,'g',lw='3')
        self.draw_line(z,y,x,R,t,'r',lw='3')

        # Inertial Frame
        R = [[1,0,0],[0,1,0],[0,0,1]]
        t = [0,0,0]
        
        self.draw_line(x,y,z,R,t,'b')
        self.draw_line(y,x,z,R,t,'g')
        self.draw_line(z,y,x,R,t,'r')

    def create_plane(self,eul = [0,0,0],t = [0,0,0]):
        ''' Create plane drawing and rotate accoring to euler angle inputs in eul
        Translate plan according to inputs in t'''
        # Get rotation matrix from euler angles
        R = self.get_rotation_matrix(eul) 

        # Clear previous drawing and label axes
        self.ax.clear() 
        self.ax.set_xlabel('X (m)')
        self.ax.set_ylabel('Y (m)')
        self.ax.set_zlabel('Z (m)')

        # Create plane
        self.create_axes(R,t)
        self.create_wing(R,t)
        self.create_fuselage(R,t)
        self.create_nose(R,t)
        self.create_tailwing(R,t)
        self.create_tailfin(R,t)

        # Draw plane
        # plt.draw()
        self.fig.canvas.draw()
        
    def get_rotation_matrix(self,eul):
        psi = eul[2] # yaw 
        th = eul[1] # pitch 
        phi = eul[0] # roll
        R = [[np.cos(th) * np.cos(psi),np.sin(phi) * np.sin(th) * np.cos(psi) - np.cos(phi) * np.sin(psi),np.cos(phi) * np.sin(th) * np.cos(psi) + np.sin(phi) * np.sin(psi)],
             [np.cos(th) * np.sin(psi),np.sin(phi) * np.sin(th) * np.sin(psi) + np.cos(phi) * np.cos(psi),np.cos(phi) * np.sin(th) * np.sin(psi) - np.sin(phi) * np.cos(psi)],
             [-np.sin(th), np.sin(phi) * np.cos(th), np.cos(phi) * np.cos(th)]]
        return R

    #--------------------------------------------------------------------------#
    # These functions do the animation post calculation
    def animate(self,i):
        ''' animation function'''
        if i < len(self.ang[0]):
            eul = [self.ang[0][i],self.ang[1][i],self.ang[2][i]]
            t = [self.pos[0][i],self.pos[1][i],self.pos[2][i]]
            self.create_plane(eul,t)

            # Edits Axes So I don't have to zoom out every time. Not sure if I want it
            if abs(t[2]) > self.bounds+5 or abs(t[1]) > self.bounds+5 or abs(t[0]) > self.bounds+5: 
                self.bounds += 10
                self.ax.set_xbound(-self.bounds,self.bounds)
                self.ax.set_ybound(-self.bounds,self.bounds)
                self.ax.set_zbound(-self.bounds,self.bounds)
            
    def show_animation(self,n_frames = 100,rep = False,interv = 10):
        self.make_figure(False)
        ani = animation.FuncAnimation(self.fig,self.animate,self.n_frames,interval=self.interv,repeat=rep)
        plt.show(block=False)
        plt.close()

    #-------------------------------------------------------------------------------#
    # These functions are for drawing real time
    def draw_update(self,eul,pos):
        '''Input euler angles and positions one time step at a time to plot real-time'''
        # Edits Axes To follow Plane around at close range
        self.ax.set_xbound(-self.bound_list[0] + pos[0],self.bound_list[0] + pos[0])
        self.ax.set_ybound(-self.bound_list[1] + pos[1],self.bound_list[1] + pos[1])
        self.ax.set_zbound(-self.bound_list[2] + pos[2],self.bound_list[2] + pos[2])
        
        # Plot
        self.create_plane(eul,pos)
        
    # # ------------------------------------------------------------------------------#
    # # If inputting post process calculated euler angles and position
    # def input_goal(self,euler,positions,single_position = False):
    #     '''Inputs:
    #     euler_angles: vector of euler angles in order psi, theta, phi
    #     positions: vector of positions in order x, y, z'''
    #     self.phi = euler[0]
    #     self.th = euler[1]
    #     self.psi = euler[2]
    #     self.x = positions[0]
    #     self.y = positions[1]
    #     self.z = positions[2]

    #     if single_position == True:
    #         self.slider_goal()
    #     else:
    #         self.calculated_trajectory()

    # def calculated_trajectory(self):
    #     self.ang = [self.phi,self.th,self.psi]
    #     self.pos = [self.x,self.y,self.z]
    #     self.show_animation()
    #     # self.plot_update()

    #     # Now set init as current
    #     self.init_psi = deepcopy(self.psi)
    #     self.init_th = deepcopy(self.th)
    #     self.init_phi = deepcopy(self.phi)
    #     self.init_x = deepcopy(self.x)
    #     self.init_y = deepcopy(self.y)
    #     self.init_z = deepcopy(self.z)


    # def slider_goal(self):
    #     # Desired euler angles and positions here (from sliders)
    #     self.ang = [np.linspace(self.init_phi,self.phi,self.n_frames),
    #                 np.linspace(self.init_th,self.th,self.n_frames),
    #                 np.linspace(self.init_psi,self.psi,self.n_frames)]
    #     self.pos = [np.linspace(self.init_x,self.x,self.n_frames),
    #                 np.linspace(self.init_y,self.y,self.n_frames),
    #                 np.linspace(self.init_z,self.z,self.n_frames)]

    #     # Run Animation
    #     self.show_animation()
                    
    #     # Now set init as current
    #     self.init_psi = deepcopy(self.psi)
    #     self.init_th = deepcopy(self.th)
    #     self.init_phi = deepcopy(self.phi)
    #     self.init_x = deepcopy(self.x)
    #     self.init_y = deepcopy(self.y)
    #     self.init_z = deepcopy(self.z)
                    

# # Create aircraft Instance
# A = AirCraftDrawing()

# # Loop if using Tkinter
# A.sm.mainloop()
    

    
