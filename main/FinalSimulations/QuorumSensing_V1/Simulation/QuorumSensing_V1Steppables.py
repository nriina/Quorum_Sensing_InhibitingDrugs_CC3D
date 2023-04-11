from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup

from cc3d.core.PySteppables import *
import numpy as np




class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)
        
        
    def start(self):
        
        for cell in self.cell_list:

            cell.targetVolume = 25
            cell.lambdaVolume = 2.0
            cell.targetSurface = 25
            cell.lambdaSurface = 2.0
            cell.dict['sensor'] = 's'
            cell.dict['sensornum'] = 0
            cell.dict['ligand'] = 'l'
            cell.dict['Me'] = 10
        
        
class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)
        self.track_cell_level_scalar_attribute(field_name='Me', attribute_name='Me')
        self.track_cell_level_scalar_attribute(field_name='sensor', attribute_name='sensornum')
        
    
    def start(self):
        
        self.plot_win = self.add_new_plot_window(title='Target Volume',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_win.add_plot("Vol", style='Dots', color='red', size=5)
        self.plot_win2 = self.add_new_plot_window(title='me',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_win2.add_plot("me", style='Dots', color='red', size=5)


    def step(self, mcs):
        
        decay_rate = 1
        mutation_prob = 0.5
        
        for cell in self.cell_list:
            s_count = 0
            if cell.dict['Me'] > 0:
                for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                        if neighbor:
                            if neighbor.dict['sensor'] == 'S':
                                s_count +=1

                if s_count > 0:
                    decay = decay_rate / s_count
                else:
                    decay = decay_rate
                cell.dict['Me'] = cell.dict['Me'] - decay
            self.plot_win2.add_data_point("me", mcs, cell.dict['Me'])

        
        if not mcs%10:
            for cell in self.cell_list:
                cell.targetVolume += 1 
                
                if np.random.uniform(0.0,1.0) < mutation_prob:
                    if cell.dict['sensor'] == 's':
                        cell.dict['sensor'] = 'S'
                        cell.dict['sensornum'] = 1
                    else:
                        cell.dict['sensor'] = 's'
                        cell.dict['sensornum'] = 0
                
                self.plot_win.add_data_point("Vol", mcs, cell.targetVolume)
                
                
                

        # # alternatively if you want to make growth a function of chemical concentration uncomment lines below and comment lines above        

        # field = self.field.CHEMICAL_FIELD_NAME
        
        # for cell in self.cell_list:
            # concentrationAtCOM = field[int(cell.xCOM), int(cell.yCOM), int(cell.zCOM)]

            # # you can use here any fcn of concentrationAtCOM
            # cell.targetVolume += 0.01 * concentrationAtCOM       

        
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list:
            if cell.volume>50:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # reducing parent target volume
        self.parent_cell.targetVolume /= 2.0                  

        self.clone_parent_2_child()            

        # for more control of what gets copied from parent to child use cloneAttributes function
        # self.clone_attributes(source_cell=self.parent_cell, target_cell=self.child_cell, no_clone_key_dict_list=[attrib1, attrib2]) 
        
        if self.parent_cell.type==1:
            self.child_cell.type=2
        else:
            self.child_cell.type=1

        
class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        for cell in self.cell_list:
            if cell.dict['Me'] <= 0:
                cell.targetVolume=0
                cell.lambdaVolume=1000

        