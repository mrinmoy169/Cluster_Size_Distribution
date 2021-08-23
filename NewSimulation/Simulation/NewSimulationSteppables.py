from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import random
import numpy as np

class NewSimulationSteppable(SteppableBasePy):

    def __init__(self,frequency=1):

        SteppableBasePy.__init__(self,frequency)
        self.force_mod = -40.01  #value of mu

    def start(self):
        list_EMT=['E','EM','M']
        probability=[0.6,0.4,0.0] #60% E cells and 40% EM (or H) cells
        k=1
        for cell in self.cell_list:
            k+=1
        def biase(lst,probability):
            zipped = zip(lst,probability)
            lst = [[i[0]] * int(i[1]*k) for i in zipped]
            new = [b for i in lst for b in i]
            return new
        ini_list = biase(list_EMT,probability)
        random.shuffle(ini_list)
        
        for cell in self.cell_list:
            if ini_list[cell.id-1] == 'E':
                cell.type = self.E
            elif ini_list[cell.id-1] == 'EM':
                cell.type = self.EM
            else:
                cell.type = self.M

            if cell.type == self.E:
                cell.targetVolume = 25
#                 cell.targetSurface = 25
            elif cell.type == self.EM:
                cell.targetVolume = 25
#                 cell.targetSurface = 25
            else:
                cell.targetVolume = 25
#                 cell.targetSurface = 25
            cell.lambdaVolume = 1.0
#             cell.lambdaSurface = 1.0
            
            #Initialize the random polarity to the cells
            if cell.type == self.EM or cell.type == self.M:
                angle = np.random.normal(0,2*np.pi)
                cell.lambdaVecX= self.force_mod*np.cos(angle)
                cell.lambdaVecY= self.force_mod*np.sin(angle)

            cell.dict["postauX"]=cell.xCOM
            cell.dict["postauY"]=cell.yCOM
            cell.dict["vx"]=0.0
            cell.dict["vy"]=0.0
     
    def step(self,mcs):
        for cell in self.cell_list:

            if mcs%50==0:    # value of tau
                dx=(cell.xCOM - cell.dict["postauX"])
                dy=(cell.yCOM - cell.dict["postauY"])
                if dx>0.5*self.dim.x:
                    dx-=self.dim.x
                if dx<=-0.5*self.dim.x:
                    dx+=self.dim.x
                if dy>0.5*self.dim.y:
                    dy-=self.dim.y
                if dy<=-0.5*self.dim.y:
                    dy+=self.dim.y
                vx = dx/50
                vy = dy/50
                cell.dict["postauX"]=cell.xCOM
                cell.dict["postauY"]=cell.yCOM
                cell.dict["vx"]=vx
                cell.dict["vy"]=vy
                norm_v = np.linalg.norm([vx,vy])
                if norm_v != 0:
                    vx /= norm_v
                    vy /= norm_v

                if cell.type==self.EM or cell.type == self.M:
                    cell.lambdaVecX = self.force_mod*vx
                    cell.lambdaVecY = self.force_mod*vy

#calculation of no. of clusters
        cluster1=[]
        cluster2=[]
        cluster3=[]
        for cell in self.cell_list_by_type(self.EM):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor is not None and neighbor.type == self.EM: 
                    cl=([cell.id,neighbor.id])
                    cluster1.append(cl)
                if neighbor is None:
                    cluster2.append(cell.id)
                if neighbor is not None and neighbor.type == self.E:
                    cluster3.append(cell.id)

        taken=[False]*len(cluster1)
        cluster_1=[set(elem) for elem in cluster1]

        def dfs(node,index):
            taken[index]=True
            ret=node
            for i,item in enumerate(cluster_1):
                if not taken[i] and not ret.isdisjoint(item):
                    ret.update(dfs(item,i))
            return ret

        ret=[]
        for i,node in enumerate(cluster_1):
            if not taken[i]:
                ret.append(list(dfs(node,i)))
        
        cluster1_separated=[item for item in ret if all(x not in cluster3 for x in item)]
        
        ret_flat=[item for sublist in ret for item in sublist]+cluster3
        cluster2_unique=[[x] for x in cluster2 if x not in ret_flat]

        cluster_EM=[len(x) for x in cluster1_separated]+[len(x) for x in cluster2_unique]
        
        #saving the data to files
        if mcs%10==0:
            output_dir = self.output_dir
            if output_dir is not None:
                output_path = Path(output_dir).joinpath('cluster_' + str(mcs) + '.dat')
                with open(output_path, 'w') as fout:
                    for item in cluster_EM:
                        fout.write('{}\n'.format(item))
        
        if mcs%10==0:
            output_dir = self.output_dir
            if output_dir is not None:
                output_path = Path(output_dir).joinpath('step_' + str(mcs) + '.dat')
                with open(output_path, 'w') as fout:
                    for cell in self.cell_list:
                        fout.write('{} {} {} {} {}\n'.format(cell.id, cell.type, cell.xCOM, cell.yCOM, cell.volume))
        

         