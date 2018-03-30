import time
from pathManager import *

class pathPlanner(pathManager):
    def __init__(self,va0 = 30.0,model=1):
        pathManager.__init__(self,va0,model)

        self.map_width = 5000
        self.map_height = 5000
        self.map_space = 500
        self.map_x = np.arange(self.map_space,self.map_width,self.map_space)
        self.map_lx = np.arange(self.map_space,self.map_width,self.map_space)
        self.map_ly = np.arange(self.map_space,self.map_width,self.map_space)
        self.map_x = np.tile(self.map_x,(len(self.map_x),1))
        self.map_y = np.transpose(self.map_x)
        self.building_width = 20
        self.map_heights = np.ones(np.shape(self.map_y)) * np.random.randn() * 30

    def planRRT(self,ps,pe):
        segmentLength = 100
        chi = 0
        # N E D chi cost parent idx flag_connect_to_goal
        start_node = np.reshape([ps[0], ps[1], pe[2], chi, 0, 0, 0],(1,7))
        end_node = np.reshape([pe[0], pe[1], pe[2], chi, 0, 0, 0],(1,7))

        tree = start_node
        if np.linalg.norm(np.subtract(start_node[0:3],end_node[0:3])) < segmentLength and self.collision(start_node,end_node) == 0:
            path = [start_node, end_node];
        else:
            numPaths = 0
            while numPaths < 2:
                tree, flag = self.extendTree(tree,end_node,segmentLength,pe[2],chi)
                numPaths = numPaths + flag
                print numPaths

        path = self.findMinimumPath(tree,end_node[0])
        print path
        path_out = self.smoothPath(path)
        print path_out
        # PLOT?

    def generateRandomNode(self,pd,chi):
        pn = self.map_width * abs(np.random.randn(1))
        pe = self.map_width * abs(np.random.randn(1))
        cost = 0
        return [pn.tolist()[0], pe.tolist()[0], pd, chi, cost, 0, 0]

    def collision(self,start_node,end_node):
        collision_flag = 0
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
                connectingNodes.append(tree[i,:])
        tmp = min(connectingNodes[:,4])
        idx = np.argmin(connectingNodes[:,4])
        
        path = [connectingNodes[idx,:], end_node]
        parent_node = connectingNodes[idx,5]
        
        while parent_node > 1:
            parent_node = tree[parent_node,5]
            path = [tree[parent_node,:], path]
            
        return path

    def smoothPath(self,path):
        newPath = path[0,:]
        ptr = 2
        while ptr <= np.size(path,1) - 1:
            if self.collision(newPath[-1,:], path[ptr+1,:]) != 0:
                newPath = [newPath, path[ptr,:]]
            ptr += 1
        newPath = [newPath, path[-1,:]]
        return newPath
