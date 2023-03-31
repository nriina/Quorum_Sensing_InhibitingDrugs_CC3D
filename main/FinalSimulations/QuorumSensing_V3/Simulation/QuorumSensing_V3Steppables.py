
from cc3d.core.PySteppables import *


## parameters
l_decay = False
s_decay = True
cheating = True

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):

        for cell in self.cell_list:

            cell.targetVolume = 40
            cell.lambdaVolume = 2.0
            cell.targetSurface = 25
            cell.lambdaSurface = 2.0
            if np.random.uniform(0.0,1.0) < 0.2:
                cell.dict['sensor'] = 'S'
                cell.dict['sensornum'] = 2
            else:
                cell.dict['sensor'] = 'S'
                cell.dict['sensornum'] = 0
            cell.dict['ligand'] = 'l'
            cell.dict['Me'] = 10
            if cheating:
                if np.random.uniform(0.0,1.0) < 0.1:
                    cell.dict['cheat'] = 1
                else:
                    cell.dict['cheat'] = 0
            else:
                cell.dict['cheat'] = 0
        
        
class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)
        
    def start(self):
        self.track_cell_level_scalar_attribute(field_name='Me', attribute_name='Me')
        self.track_cell_level_scalar_attribute(field_name='sensor', attribute_name='sensornum')
        self.track_cell_level_scalar_attribute(field_name='cheating', attribute_name='cheat')
        
        self.plot_win = self.add_new_plot_window(title='Target Volume',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_win.add_plot("Vol", style='Dots', color='red', size=5)
        self.plot_win2 = self.add_new_plot_window(title='Metabolic Energy',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Volume', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_win2.add_plot("me", style='Dots', color='red', size=5)
        
        # Make sure Secretion plugin is loaded
        # make sure this field is defined in one of the PDE solvers
        # you may reuse secretor for many cells. Simply define it outside the loop
        
        
        # # arguments are: cell, max uptake, relative uptake
        # secretor.uptakeInsideCellAtBoundary(cell, 2.0, 0.2)
        # # tot_amount = secretor.uptakeInsideCellAtBoundaryTotalCount(cell, 2.0, 0.2).tot_amount
        
        

    def step(self, mcs):
    
        #this code does metabolic energy with a decay, 
        #which is modulated by the amount of neighbors that are s allele
        
        secretor_ai = self.get_field_secretor("ai")
        secretor_product = self.get_field_secretor("product")
        secretor_oxy = self.get_field_secretor("oxy")
        
        decay_rate = 0.5
        mutation_prob = 0.05
        ai_thresh = 30
        
        for cell in self.cell_list:
            tot_amount_ai = secretor_ai.uptakeInsideCellAtBoundaryTotalCount(cell, 2.0, 0.2).tot_amount
            tot_amount_product = secretor_product.uptakeInsideCellAtBoundaryTotalCount(cell, 2.0, 0.2).tot_amount
            
            '''
            note to nick, add linear interpolation function for product-oxygen uptake
            
            
            oxgyen max value = 0.5 = 10 oxygen
            '''
            # if product update is = 2(max):
                # update 1.0 oxygen
              # if product uptake = 0:
                  # uptake 0.2 oxygen
            
            
            oxy_max = 0.5 # = 10 oxygen at max
            # oxy_up = -tot_amount_product + 0.2
            oxy_up = 0.2
            tot_amount_oxygen = secretor_oxy.uptakeInsideCellAtBoundaryTotalCount(cell, oxy_max, oxy_up).tot_amount
            self.plot_win.add_data_point("Vol", mcs, -tot_amount_oxygen)  
            # s_count = 0
            decay = 0
            
            #metabolic energy - decay () - (constant if L) and (constant if S) + oxygen
            if cell.dict['Me'] > 0:
                if cell.dict['Me'] > 100:
                    cell.dict['Me'] = 100
                #code for neighbor influnce
                # for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                        # if neighbor:
                            # if neighbor.dict['sensor'] == 'S':
                                # s_count +=1
                                

                # if s_count > 0:
                    # decay = decay_rate / s_count
                # else:
                    # decay = decay_rate
                if s_decay:
                    if cell.dict['sensor'] == 'S':
                        decay +=1
        
                # if l_decay:    
                    # if cell.dict['ligand'] == 'L':
                        # decay +=1

                cell.dict['Me'] = cell.dict['Me'] - tot_amount_oxygen - (decay + decay_rate) 
            self.plot_win2.add_data_point("me", mcs, cell.dict['Me'])

        #this code mutates cells between s and S alleles
        if not mcs%10:
            for cell in self.cell_list:
                cell.targetVolume += 1 
                
                # Make sure Secretion plugin is loaded
                # make sure this field is defined in one of the PDE solvers
                # you may reuse secretor for many cells. Simply define it outside the loop

                secretor_ai.secreteOutsideCellAtBoundary(cell, 50)
                if cell.dict['sensor'] == 'S':
                    if not cell.dict['cheat']:
                        secretor_product.secreteOutsideCellAtBoundary(cell, 20)
                    
                if -tot_amount_ai > ai_thresh:
                    if cell.dict['sensor'] == 's' and cell.dict['sensornum'] == 0:
                        cell.dict['sensor'] = 'S'
                        cell.dict['sensornum'] = 1
                # tot_amount = secretor.secreteOutsideCellAtBoundaryTotalCount(cell, 300).tot_amount
                
                
                
                
                if np.random.uniform(0.0,1.0) < mutation_prob:
                    if cell.dict['sensor'] == 's':
                        cell.dict['sensor'] = 'S'
                        cell.dict['sensornum'] = 2
                    else:
                        cell.dict['sensor'] = 's'
                        cell.dict['sensornum'] = 0
                
                # self.plot_win.add_data_point("Vol", mcs, cell.targetVolume)     

        
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
        
        # if self.parent_cell.type==1:
            # self.child_cell.type=2
        # else:
            # self.child_cell.type=1

        
class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        for cell in self.cell_list:
            if cell.dict['Me'] <= 0:
                cell.targetVolume=0
                cell.lambdaVolume=1000

        