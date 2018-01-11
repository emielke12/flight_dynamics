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
        self.fig = plt.figure(figsize=(10,10)) # Figure Window size
        self.ax = a3(self.fig) # 3d axis
        self.ax.set_xbound(-5,5) # Axis bounds
        self.ax.set_ybound(-5,5)
        self.ax.set_zbound(-5,5)
        self.ax.view_init(elev=21,azim=43) # Axis view
        self.ax.dist=10 
        self.ax.invert_zaxis() # North east down view
        self.ax.invert_yaxis() # North east down view
        
    def draw_2dpolygon(self,px,py,pz,R,t,colour='r'):
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
        self.create_wing(R,t)
        self.create_fuselage(R,t)
        self.create_nose(R,t)
        self.create_tailwing(R,t)
        self.create_tailfin(R,t)

        # Draw plane
        plt.draw()
        
    def get_rotation_matrix(self,eul):
        psi = -eul[0] # yaw (is negative due to matplotlib limitations)
        th = -eul[1] # pitch (is negative due to matplotlib limitations)
        phi = eul[2] # roll
        R = [[np.cos(th)*np.cos(psi), np.cos(th)*np.sin(psi), -np.sin(th)],
             [np.sin(phi)*np.sin(th)*np.cos(psi) - np.cos(phi)*np.sin(psi), np.sin(phi)*np.sin(th)*np.sin(psi) + np.cos(phi)*np.cos(psi), np.sin(phi)*np.cos(th)],
             [np.cos(phi)*np.sin(th)*np.cos(psi) + np.sin(phi)*np.sin(psi), np.cos(phi)*np.sin(th)*np.sin(psi) - np.sin(phi)*np.cos(psi), np.cos(phi)*np.cos(th)]] # rotation matrix
        return R

    def animate(self,i):
        ''' animation function'''
        eul = [0,0,ang[i]]
        t = [0,0,i/20.0]
        self.create_plane(eul,t)


A = AirCraftDrawing()
ang = np.linspace(0,np.pi/2,100)
for j in xrange(90):
    ani = animation.FuncAnimation(A.fig,A.animate,len(ang),interval=10,repeat = True)
plt.show()
