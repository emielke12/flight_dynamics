import time
from pathManager import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

class pathPlanner(pathManager):
    def __init__(self,va0 = 30.0,model=1):
        pathManager.__init__(self,va0,model)

        self.map_width = 1000
        self.map_height = 1000
        self.map_space = 150
        self.map_x = np.arange(self.map_space,self.map_width,self.map_space)
        self.map_lx = np.arange(self.map_space,self.map_width,self.map_space)
        self.map_ly = np.arange(self.map_space,self.map_width,self.map_space)
        self.map_x = np.tile(self.map_x,(len(self.map_x),1))
        self.map_y = np.transpose(self.map_x)
        self.building_width = 20
        s = np.shape(self.map_y)
        self.map_heights = np.abs(np.random.randn(s[0],s[1])) * 100

    def planRRT(self,ps,pe,chi):
        segmentLength = 3 * self.Rmin
        # N E D chi cost parent idx flag_connect_to_goal
        start_node = np.reshape([ps[0], ps[1], pe[2], chi, 0, 0, 0],(1,7))
        end_node = np.reshape([pe[0], pe[1], pe[2], chi, 0, 0, 0],(1,7))

        tree = start_node
        if np.linalg.norm(np.subtract(start_node[0:3],end_node[0:3])) < segmentLength and self.collision(start_node,end_node) == 0:
            print 'Only One length'
            tree = [start_node, end_node]
            start_node = np.reshape(start_node,(7,))
            end_node = np.reshape(end_node,(7,))
            P = np.array([[start_node[0],start_node[1],start_node[2],chi],[end_node[0],end_node[1],end_node[2],0]])
        else:
            numPaths = 0
            while numPaths < 1:
                tree, flag = self.extendTree(tree,end_node,segmentLength,pe[2],chi)
                numPaths = numPaths + flag

            path = self.findMinimumPath(tree,end_node[0])
            path_out = self.smoothPath(path)
            self.plot_paths(path,path_out)
            P = path_out[:,0:4]
            ps = [ps[0],ps[1],ps[2],chi]
            P = np.concatenate((np.reshape(ps,(1,4)),P))
            for i in xrange(len(P) - 1):
                P[i+1,3] = self.getAngle(P[i,0:3],P[i+1,0:3],P[i+1,3])

        return P.tolist()

    def getAngle(self,ps,pe,chi_next):
        x = pe[0] - ps[0]
        y = pe[1] - ps[1]
        return np.arctan2(y,x) + (chi_next - np.arctan2(y,x))/2

    def generateRandomNode(self,pd,chi):
        pn = self.map_width * abs(np.random.randn(1))
        pe = self.map_width * abs(np.random.randn(1))
        cost = 0
        return [pn.tolist()[0], pe.tolist()[0], pd, chi, cost, 0, 0]

    def collision(self,start_node,end_node):
        collision_flag = 0
        start_node = np.reshape(start_node,(7,))
        end_node = np.reshape(end_node,(7,))
        X,Y,Z = self.pointsAlongPath(start_node,end_node,0.1)
        for i in xrange(len(X)):
            if Z[i] >= self.downAtNE(X[i],Y[i]):
                collision_flag = 1
        return collision_flag

    def pointsAlongPath(self,start_node,end_node,Del):
        X = [start_node[0]]
        Y = [start_node[1]]
        Z = [start_node[2]]
        q = np.subtract(end_node[0:3],start_node[0:3])
        L = np.linalg.norm(q)
        q /= L
        w = start_node[0:3]

        for i in xrange(int(np.floor(L/Del)) - 1):
            w += Del * q
            X.append(w[0])
            Y.append(w[1])
            Z.append(w[2])

        return X,Y,Z

    def downAtNE(self,n,e):
        d_n = np.min(abs(np.subtract(n,self.map_lx)))
        d_e = np.min(abs(np.subtract(e,self.map_ly)))
        idx_n = np.argmin(abs(np.subtract(n,self.map_lx)))
        idx_e = np.argmin(abs(np.subtract(e,self.map_ly)))

        if d_n <= self.building_width and d_e <= self.building_width:
            down = -self.map_heights[idx_e,idx_n]
        else:
            down = 0

        return down

    def extendTree(self,tree,end_node,segmentLength,pd,chi):
        flag1 = 0
        while flag1 == 0:
            randomNode = self.generateRandomNode(pd,chi)
            tmp = np.subtract(np.transpose(tree)[:][0:3],
                              np.transpose(np.ones((np.size(tree,0),1)) * randomNode[0:3]))
            tmp = np.transpose(tmp)
            dist = min(np.diag(np.matmul(tmp,np.transpose(tmp))))
            idx = np.argmin(np.diag(np.matmul(tmp,np.transpose(tmp))))
            L = min([np.sqrt(dist),segmentLength])
            cost = tree[idx][4] + L
            tmp = np.subtract(randomNode[0:3],tree[idx][0:3])
            new_point = np.add(tree[idx][0:3],L*(tmp / np.linalg.norm(tmp))).tolist()
            new_node = [new_point[0], new_point[1], new_point[2], chi, cost, idx, 0]

            if self.collision(tree[idx][:], new_node) == 0:
                new_tree = np.concatenate((tree,np.reshape(new_node,(1,7))))
                flag1 = 1

        if np.linalg.norm(np.subtract(new_node[0:3],end_node[0][0:3])) < segmentLength and self.collision(new_node,end_node[0]) == 0:
            flag = 1
            new_tree[-1,6] = 1
        else:
            flag = 0

        return new_tree, flag

    def findMinimumPath(self,tree,end_node):
        connectingNodes = []
        for i in xrange(np.size(tree,0)):
            if tree[i,-1] == 1:
                connectingNodes.append(tree[i,:].tolist())
        connectingNodes = np.transpose(connectingNodes)
        tmp = min(connectingNodes[4,:])
        idx = np.argmin(connectingNodes[4,:])
        connectingNodes = np.transpose(connectingNodes)
        path = [connectingNodes[idx,:].tolist(), end_node.tolist()]
        parent_node = int(connectingNodes[idx,5])
        while parent_node > 1:
            parent_node = int(tree[parent_node,5])
            path = np.concatenate((np.reshape(tree[parent_node,:],(1,7)),path))
            
        return path

    def smoothPath(self,path):
        l = np.size(path,1)
#         newPath = np.reshape(path[0,:],(1,l))
        newPath = np.reshape(path[0][:],(1,l))
        ptr = 2
        while ptr <= np.size(path,0) - 1:
            if self.collision(deepcopy(newPath[-1][:]), path[ptr][:]) != 0:
                newPath = np.concatenate((newPath,np.reshape(path[ptr-1][:],(1,l))))
            ptr += 1
        newPath = np.concatenate((newPath, np.reshape(path[-1][:],(1,l))))
        return newPath

    def plot_paths(self,path,smooth):
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(self.map_x,self.map_y,self.map_heights)
        redo_path = np.transpose(path)
        redo_smooth =  np.transpose(smooth)

#         path_x = path[:][0]
#         path_y = path[:][1]
#         path_z = -path[:][2]
#         smooth_x = smooth[:][0]
#         smooth_y = smooth[:][1]
#         smooth_z = -smooth[:][2]

        path_x = redo_path[0][:]
        path_y = redo_path[1][:]
        path_z = -redo_path[2][:]
        smooth_x = redo_smooth[0][:]
        smooth_y = redo_smooth[1][:]
        smooth_z = -redo_smooth[2][:]

        ax.plot(path_x,path_y,path_z,'--')
        ax.plot(smooth_x,smooth_y,smooth_z)
        plt.show()

    def planCover(self,p,ps,chi):
        pd = p[2]

        start_node = [p[0],p[1],pd,chi,0,0]

        returnMapSize = 30
        return_map = 50 * np.ones((returnMapSize,returnMapSize)) + np.random.randn(returnMapSize,returnMapSize)
#         self.plotReturnMap(return_map,returnMapSize)

        search_cycles = 200
        L = 50
        vartheta = np.pi/6
        depth = 5

        path = start_node
        path = np.reshape(path,(1,len(path)))
        for i in xrange(search_cycles):
            tree = self.extendTreeCover(path[-1,:],L,vartheta,depth,return_map,pd)
            next_path = self.findMaxReturnPath(tree)
            path = np.concatenate((path,np.reshape(next_path[0,:],(1,6))))
            return_map = self.updateReturnMap(np.reshape(next_path[0,:],(1,6)),return_map)
#             self.plotReturnMap(return_map)
        self.plotReturnMap(return_map,returnMapSize)

        path_ = path
        path = np.reshape(path[0,:],(1,6))
        for i in xrange(np.size(path_,0)-1):
            if path_[i+1,3] != path_[i,3]:
                path = np.concatenate((path,np.reshape(path_[i+1,:],(1,6))))

        path = self.smoothPath(path)
        P = path[:,0:4]
        ps = [ps[0],ps[1],ps[2],chi]
        for i in xrange(len(P) - 1):
            P[i+1,3] = self.getAngle(P[i,0:3],P[i+1,0:3],P[i+1,3])
        return P.tolist()

    def extendTreeCover(self,start_node,L,vartheta,depth,return_map,pd):
        tree_ = np.reshape(np.append(start_node,0),(1,7))
        for d in xrange(depth):
            newnodes = []
            for j in xrange(np.size(tree_,0)):
                if tree_[j,6] != 1:
                    for i in xrange(3):
                        if i == 0:
                            theta = tree_[j,3] - vartheta
                        elif i == 1:
                            theta = tree_[j,3]
                        elif i == 2:
                            theta = tree_[j,3] + vartheta
                        newnode_ = [tree_[j,0] + L * np.cos(theta),
                                    tree_[j,1] + L * np.sin(theta),
                                    tree_[j,2],
                                    theta,
                                    0,
                                    j,
                                    0]
                        if self.collisionCover(tree_[j,:], newnode_) == 0:
                            newnode_[4] = tree_[j,4] + self.findReturn(newnode_[0],newnode_[1],return_map)
                            newnodes = np.concatenate((newnodes,newnode_))
                    tree_[j,6] = 1
            tree_ = np.concatenate((tree_,np.reshape(newnodes,(-1,7))))
        tree = tree_[:,:6]
        return tree

    def findMaxReturnPath(self,tree):
        tmp = np.max(tree[:,4])
        idx = np.argmax(tree[:,4])

        path = np.reshape(tree[idx,:],(1,6))
        parent_node = tree[idx,5]
        while parent_node > 1:
            path = np.concatenate((np.reshape(tree[int(parent_node),:],(1,6)),path))
            parent_node = tree[int(parent_node),5]
        return path

    def collisionCover(self,start_node,end_node):
        collision_flag = 0
        sigma = np.linspace(0,1,11)
        for i in xrange(11):
            X = (1 - sigma[i]) * start_node[0] + sigma[i] * end_node[0]
            Y = (1 - sigma[i]) * start_node[1] + sigma[i] * end_node[1]
            Z = (1 - sigma[i]) * start_node[2] + sigma[i] * end_node[2]
            if Z >= self.downAtNE(X,Y):
                collision_flag = 1
                
            if X > self.map_width or X < 0 or Y > self.map_width or Y < 0:
                collision_flag = 1

        return collision_flag

    def findReturn(self,pn,pe,return_map):
        pnsize = np.shape(return_map)
        pn_max = pnsize[0]-1
        pe_max = pnsize[1]-1

        fn = pn_max * pn / self.map_width
        fn = min(pn_max,int(round(fn)))
        fn = max(1,fn)
        fe = pe_max * pe / self.map_width
        fe = min(pe_max,int(round(fe)))
        fe = max(1,fe)
        return return_map[fn,fe]
        

    def updateReturnMap(self,path,return_map):
        new_return_map = return_map
        for i in xrange(np.size(path,0)):
            pn = path[i,0]
            pe = path[i,1]
            pnsize = np.shape(return_map)
            pn_max = pnsize[0]-1
            pe_max = pnsize[1]-1
            fn = pn_max * pn / self.map_width
            fn = min(pn_max, int(round(fn)))
            fn = max(1,fn)
            fe = pe_max * pe / self.map_width
            fe = min(pe_max, int(round(fe)))
            fe = max(1,fe)
            
            new_return_map[fn,fe] = return_map[fn,fe] - 50

        return new_return_map


    def plotReturnMap(self,return_map,s):
        self.fig_surface = plt.figure()
        self.ax_surface = self.fig_surface.add_subplot(111,projection='3d')  
        l = np.size(return_map,0)
        X = np.arange(0,l,l/s)
        Y = np.arange(0,l,l/s)
        X,Y = np.meshgrid(X,Y)
        Z = return_map
        self.ax_surface.plot_surface(X,Y,Z,cmap = cm.coolwarm)
        plt.show()
        

        
