from mpl_toolkits.mplot3d import Axes3D as a3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection as pc
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time,sys
from copy import deepcopy
from Tkinter import *

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
        self.axes_length = 1.5

        # Figure
        plt.close('all')
        self.fig = plt.figure(figsize=(10,10)) # Figure Window size
        self.ax = a3(self.fig) # 3d axis
        self.ax.set_xbound(-5,5) # Axis bounds
        self.ax.set_ybound(-5,5)
        self.ax.set_zbound(-5,5)
        self.ax.view_init(elev=21,azim=43) # Axis view
        self.ax.dist=10 
        self.ax.invert_zaxis() # North east down view
        self.ax.invert_yaxis() # North east down view

        # Animation Values
        self.n_frames = 40
        self.interv = 10

        # Loop condition
        self.loop_start = False
        self.loop = True

        # Slider Definition
        self.sm = Tk()
        self.phi_s = Scale(self.sm, from_=0, to=360,orient=HORIZONTAL,resolution=5,label='Roll')
        self.th_s = Scale(self.sm, from_=0, to=360,orient=HORIZONTAL,resolution=5,label='Pitch')
        self.psi_s = Scale(self.sm, from_=0, to=360,orient=HORIZONTAL,resolution=5,label='Yaw')
        self.x_s = Scale(self.sm, from_=-10, to=10,orient=HORIZONTAL,resolution=0.5,label='X')
        self.y_s = Scale(self.sm, from_=-10, to=10,orient=HORIZONTAL,resolution=0.5,label='Y')
        self.z_s = Scale(self.sm, from_=-10, to=10,orient=HORIZONTAL,resolution=0.5,label='Z')
        self.phi_s.pack()
        self.th_s.pack()
        self.psi_s.pack()
        self.x_s.pack()
        self.y_s.pack()
        self.z_s.pack()
        self.phi_s.set(0)
        self.th_s.set(0)
        self.psi_s.set(0)
        self.x_s.set(0)
        self.y_s.set(0)
        self.z_s.set(0)

        # Button Definition
        Button(self.sm,text='Plot',command=self.start).pack()
        Button(self.sm,text='Reset Sliders',command=self.reset).pack()
        Button(self.sm,text='Stop',command=self.stop).pack()
        
        # Visualization start points
        self.phi = 0.0
        self.th = 0.0
        self.psi = 0.0
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.init_phi = 0.0
        self.init_th = 0.0
        self.init_psi = 0.0
        self.init_x = 0.0
        self.init_y = 0.0
        self.init_z = 0.0

        
    def draw_line(self,px,py,pz,R,t,colour='r',lw = '5'):
        ''' Draw a line for axes'''
        new = np.split(np.dot(R,np.column_stack((np.column_stack((px,py)),pz)).T),3)
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1])
        z = np.add(new[2][0],t[2])
        verts = [zip(x,y,z)]
        self.ax.add_collection3d(pc(verts,edgecolors=colour,linewidth=lw))

        
    def draw_2dpolygon(self,px,py,pz,R,t,colour='#929591'):
        '''Draw a polygon with inputs px,py,pz which are a list of vertices
        Rotate the points according to the specified rotation matrix R
        Translate the points according to the translation vector t
        Give polygon a color according to 'r' '''
        # # Uncomment to see erroneous translation before rotation
        # px = np.add(px,t[0])
        # py = np.add(py,t[1])
        # pz = np.add(pz,t[2])

        new = np.split(np.dot(R,np.column_stack((np.column_stack((px,py)),pz)).T),3)

        # Uncomment to see correct rotation before translation
        x = np.add(new[0][0],t[0])
        y = np.add(new[1][0],t[1])
        z = np.add(new[2][0],t[2])
        
        # # Uncomment to see erroneous translation before rotation
        # x = new[0][0]
        # y = new[1][0]
        # z = new[2][0]

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
        R = self.get_rotation_matrix(eul) # Get rotation matrix from euler angles

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
        plt.draw()
        
    def get_rotation_matrix(self,eul):
        psi = eul[0] # yaw 
        th = eul[1] # pitch 
        phi = eul[2] # roll
        R = [[np.cos(th)*np.cos(psi), np.cos(th)*np.sin(psi), -np.sin(th)],
             [np.sin(phi)*np.sin(th)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(th)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(th)],
             [np.cos(phi)*np.sin(th)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(th)]] # rotation matrix
        R = np.transpose(R)
        # R_roll = [[1,0,0],
        #           [0,np.cos(phi),np.sin(phi)],
        #           [0,-np.sin(phi),np.cos(phi)]]
        # R_pitch = [[np.cos(th),0,-np.sin(th)],
        #           [0,1,0],
        #           [np.sin(th),0,np.cos(th)]]
        # R_yaw = [[np.cos(psi),np.sin(psi),0],
        #           [-np.sin(psi),np.cos(psi),0],
        #           [0,0,1]]
        # R = np.transpose(np.matmul(np.matmul(R_roll,R_pitch),R_yaw))
        return R

    #--------------------------------------------------------------------------#
    # These functions do the animation
    def animate(self,i):
        ''' animation function'''
        eul = [self.ang[0][i],self.ang[1][i],self.ang[2][i]]
        t = [self.pos[0][i],self.pos[1][i],self.pos[2][i]]
        self.create_plane(eul,t)

    def show_animation(self,n_frames = 100,rep = False,interv = 10):
        ani = animation.FuncAnimation(self.fig,self.animate,self.n_frames,interval=self.interv,repeat=rep)
        plt.show(block=False)
        plt.close()

    # ------------------------------------------------------------------------------#
    # These functions deal with the Tkinter block
    def get_values(self):
        # Gets values from Tkinter sliders
        self.psi = self.psi_s.get() * np.pi/180.0
        self.th = self.th_s.get() * np.pi/180.0
        self.phi = self.phi_s.get() * np.pi/180.0
        self.x = self.x_s.get()
        self.y = self.y_s.get()
        self.z = self.z_s.get()

    def start(self):
        # Starts animation on push button input
        self.get_values()
        self.slider_goal()
        self.show_animation()

        # Figure
        plt.close('all')
        self.fig = plt.figure(figsize=(10,10)) # Figure Window size
        self.ax = a3(self.fig) # 3d axis
        self.ax.set_xbound(-5,5) # Axis bounds
        self.ax.set_ybound(-5,5)
        self.ax.set_zbound(-5,5)
        self.ax.view_init(elev=21,azim=43) # Axis view
        self.ax.dist=10 
        self.ax.invert_zaxis() # North east down view
        self.ax.invert_yaxis() # North east down view


    def stop(self):
        # Stops script
        sys.exit()

    def reset(self):
        # Resets slider values and position
        self.init_psi = 0.0
        self.init_th = 0.0
        self.init_phi = 0.0
        self.init_x = 0.0
        self.init_y = 0.0
        self.init_z = 0.0
        self.phi_s.set(0)
        self.th_s.set(0)
        self.psi_s.set(0)
        self.x_s.set(0)
        self.y_s.set(0)
        self.z_s.set(0)

        
    def slider_goal(self):        
        # Desired euler angles and positions here (from sliders)
        self.ang = [np.linspace(self.init_psi,self.psi,self.n_frames),
               np.linspace(self.init_th,self.th,self.n_frames),
               np.linspace(self.init_phi,self.phi,self.n_frames)]
        self.pos = [np.linspace(self.init_x,self.x,self.n_frames),
               np.linspace(self.init_y,self.y,self.n_frames),
               np.linspace(self.init_z,self.z,self.n_frames)]

        # Now set init as current
        self.init_psi = deepcopy(self.psi)
        self.init_th = deepcopy(self.th)
        self.init_phi = deepcopy(self.phi)
        self.init_x = deepcopy(self.x)
        self.init_y = deepcopy(self.y)
        self.init_z = deepcopy(self.z)

    # ------------------------------------------------------------------------------#
    # If not using Tkinter but just inputting calculated euler angles and position
    def input_goal(self,euler_angles,positions):
        '''Inputs:
        euler_angles: vector of euler angles in order psi, theta, phi
        positions: vector of positions in order x, y, z'''
        self.ang = [euler[0],euler[1],euler[2]]
        self.pos = [positions[0],positions[1],positions[2]]
        self.show_animation()



# Create aircraft Instance
A = AirCraftDrawing()

# Loop if using Tkinter
A.sm.mainloop()
    

    
