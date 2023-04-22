from cc3d.cpp.PlayerPython import * 
from cc3d import CompuCellSetup
import numpy as np
from cc3d.core.PySteppables import *

# Simulation Parameters
mutation = True
s_decay_modifier = True
cheating_consequence = True
starvation_bonus = True
wt_invasion = True
qsi_active = True

dilution = True
distance_btw_dilutions = 1000
dilution_bottleneck = True
dilution_factor = 0.75
track_cooperators = True


class ParamSteeringSteppable(SteppableBasePy):
    def __init__(self, frequency=10):
        SteppableBasePy.__init__(self, frequency)

    def add_steering_panel(self):
        self.add_steering_param(name='mutation', val=1, enum=[1,0], widget_name='combobox')
        self.add_steering_param(name='qs_decay', val=1, enum=[1,0], widget_name='combobox')
        self.add_steering_param(name='cheating_consequence', val=1, enum=[1,0], widget_name='combobox')
        self.add_steering_param(name='starvation_bonus', val=1, enum=[1,0], widget_name='combobox')
        self.add_steering_param(name='wt_invasion', val=1, enum=[1,0], widget_name='combobox')
        self.add_steering_param(name='qsi_active', val=1, enum=[1, 0], widget_name='combobox')
        
        self.add_steering_param(name='dilution', val=1, enum=[1, 0], widget_name='combobox')
        self.add_steering_param(name='distance_btw_dilutions', val=1000)
        self.add_steering_param(name='dilution_bottleneck', val=1, enum=[1, 0], widget_name='combobox')
        self.add_steering_param(name='dilution_factor', val=0.75)
        self.add_steering_param(name='track_cooperators', val=1, enum=[1, 0], widget_name='combobox')
        self.add_steering_param(name='mutation_probability', val=0.001)

    def process_steering_panel_data(self):
        global mutation
        global s_decay_modifier
        global cheating_consequence
        global starvation_bonus
        global wt_invasion
        global qsi_active
        
        global dilution
        global distance_btw_dilutions
        global dilution_bottleneck
        global dilution_factor
        global track_cooperators
        
        mutation = self.get_steering_param('mutation')
        s_decay_modifier = self.get_steering_param('qs_decay')
        cheating_consequence = self.get_steering_param('cheating_consequence')
        starvation_bonus = self.get_steering_param('starvation_bonus')
        wt_invasion = self.get_steering_param('wt_invasion')
        qsi_active = self.get_steering_param('qsi_active')
        
        dilution = self.get_steering_param('dilution')
        distance_btw_dilutions = self.get_steering_param('distance_btw_dilutions')
        dilution_bottleneck = self.get_steering_param('dilution_bottleneck')
        dilution_factor = self.get_steering_param('dilution_factor')
        track_cooperators = self.get_steering_param('track_cooperators')
        

"""
mutation parameter decides whether or not cells will mutate as the simulation progresses. Its
probability can be modulated with the mutation_prob local variable within the ConstraintInitializer
Steppable step funciton.

s_decay_modifier parameter decides whether or not there is a metabolic consequence to active
Quorum Sensing. This adds the to the cell type specific decay modifier for any quorum sensing cells.

cheating_consequence parameter decides whether or not defector cell types get a penalty for not 
quorum sensing. This penalty is an increase in Me decay that has an exponential relationship to
the concentration of autoinducer (AI). The try except statement is because the AI conc. can cause
a value that is too large. The except just sets this to the max set within the try statement.

starvation_bonus parameter decreases the energy decay of cells that are QS capable. It also
decreases growth. This only comes into effect if the cells are starving and running out of stored
Me

wt_invasion parameter describes whether or not a wild-type cell is spawned after 500 mcs.
This is only used to test how majority UC or D populations respond to the appearance of WT

qsi_active parameter describes whether or not to include the effects of the QSI chemical field.
The IC50 of any quorum sensing inhibitor can be simulated by changing the reciprocal function
that describes how the public good secretion is effected by the drug. Secretion of 1 = 100% 
Secretion of 0.5 =  50% QS activitiy.

Dilution refers to the broad set of parameters that deals with addition of fresh media (nutrients)

distance_btw _dilution refers to the amount of mcs between dilution events

dilution_bottleneck True means that cells are randomly killed to simulate backdiltion from a larger population

dilution_factor refers to the percentage of cells that are killed in a dilution event.

track_cooperators produces a figure that outputs the fraction of population that is cooperating (QS)
"""
# Function Library
def modify_value(variable, ceiling, modifier, exponential=True, linear=False):

    ceiling = ceiling * variable
    
    if linear and not exponential:
        if modifier < ceiling:
            mod_var = variable + modifier
        else:
            mod_var = variable + ceiling
        
        return mod_var
    
    elif exponential and not linear:
        if modifier < ceiling:
            mod_var = variable * modifier
        else:
            mod_var = variable * modifier
        
        return mod_var
    
    else:
        print('Modification cannot be both exponential and linear. Set one to False')
                
class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self,frequency)

    def start(self):

        for cell in self.cell_list:

            cell.targetVolume = 40
            cell.lambdaVolume = 4.0
            
            cd = self.chemotaxisPlugin.addChemotaxisData(cell, "complex nutrient")
            cd.setLambda(30.0)
            
            nd = self.chemotaxisPlugin.addChemotaxisData(cell, 'nutrient')
            nd.setLambda(20.0)


            
            # Create Initial Sensed AI key 
            cell.dict['Sensed AI'] = 0
            
            # Create Key of Metabolic Energy 'Me' for each cell
            cell.dict['Me'] = 10
            
            # Set basal quorum sensing sensor state for each cell type 1 - on; 0 - off
            if cell.type == 1: # Wild-Type (WT) - Color is Blue
                cell.dict['sensor'] = 0
                
            elif cell.type == 2: # Unconditional Cooperator (UC) - Color is Green
                cell.dict['sensor'] = 1
                
            elif cell.type == 3: # Defector (D) - Color is Red
                cell.dict['sensor'] = 0
        
        
class GrowthSteppable(SteppableBasePy):
    def __init__(self,frequency=1):
        SteppableBasePy.__init__(self, frequency)


    def start(self):
        
        # Sensor Activity Field
        self.track_cell_level_scalar_attribute(field_name='Sensor', attribute_name='sensor') 
        
        # Metabolic Energy Tracking Field
        self.track_cell_level_scalar_attribute(field_name='Me', attribute_name='Me')
        
        # Plot Metabolic Energy
        self.plot_me = self.add_new_plot_window(title='Metabolic Energy',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_me.add_plot("WT Me", style='Dots', color='blue', size=3)
        self.plot_me.add_plot("UC Me", style='Dots', color='green', size=3)
        self.plot_me.add_plot("D Me", style='Dots', color='red', size=3)
        
        # Plot Cell Type Populations
        self.plot_cell_pops = self.add_new_plot_window(title='Fraction of Population by Cell Type',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='% of Population', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_cell_pops.add_plot('WT', style='Lines', color='blue', size=3)
        self.plot_cell_pops.add_plot('UC', style='Lines', color='green', size=3)
        self.plot_cell_pops.add_plot('D', style='Lines', color='red', size=3)
        
        # Plot AI / QSI Uptake
        self.plot_ai = self.add_new_plot_window(title='AI Uptake',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
        
        self.plot_ai.add_plot("AI Uptake", style='Lines', color='purple', size=3)
        
        if track_cooperators:
            self.plot_coop = self.add_new_plot_window(title='Frequency of Cooperation',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Frequency', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
            
            self.plot_coop.add_plot('Cooperating Fraction', style='Lines', color='white', size=3)
        
        if qsi_active:
            self.plot_qsi = self.add_new_plot_window(title='QSI Uptake',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Concentration (uM)', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False)
            
            self.plot_qsi.add_plot('QSI Uptake', style="Lines", color='yellow', size=3)
            self.plot_qsi.add_plot('Public Good Excretion', style='Lines', color='lime', size=3)

    
    def step(self, mcs):
        print("step qsi value:",qsi_active)
        
        if dilution:
            if not (mcs % distance_btw_dilutions):
                self.field.nutrient[0:100, 0:100, 0] = 20
                self.field.complex_nutrient[0:100, 0:100, 0] = 50
                self.field.ai[0:100, 0:100, 0] = 0
                self.field.public_good[0:100, 0:100, 0] = 0
        
        # WT invasion
        if wt_invasion:
            # Spawn WT After Simulation initiation 
            if mcs == 1001 and len(self.cell_list_by_type(self.D)) < 1:
                # size of cell will be 3x3x1
                self.cell_field[50:54, 50:54, 0] = self.new_cell(self.D)
                for cell in self.cell_list_by_type(self.D):
                    cell.dict['Me'] = 1000
                    cell.targetVolume = 40
                    cell.lambdaVolume = 4.0
                    cell.dict['sensor'] = 0
            
            

        # Establish Secretors
        secretor_ai = self.get_field_secretor("ai")
        secretor_public_good = self.get_field_secretor("public_good")
        secretor_nutrient = self.get_field_secretor("nutrient")
        secretor_complex = self.get_field_secretor('complex_nutrient')
        secretor_qsi = self.get_field_secretor('qsi')
        
        # Define Average AI Uptake
        raw_ai_uptake= 0
        raw_qsi_uptake = 0
        raw_public_good_excretion = 0
        
        
        # Define Function Variables
        decay_rate = 0.05
        mutation_prob = self.get_steering_param('mutation_probability')
        starvation_mutation_prob = 0.001
        ai_threshold = 0.7
        ai_lower_thres = 0.5
        nutrient_rel_uptake = 0.02
        nutrient_max_uptake = 0.5
    
        for cell in self.cell_list:
            # Field Values for each cell
            total_ai = -1*secretor_ai.uptakeInsideCellAtBoundaryTotalCount(cell, 0.5, 0.02).tot_amount 
            if qsi_active:
                total_qsi = -1 * secretor_qsi.uptakeInsideCellAtBoundaryTotalCount(cell, 0.1, 0.02 ).tot_amount
                
                if total_qsi < 0: # sets floor of qsi uptake
                    total_qsi = 0
                
                raw_qsi_uptake += total_qsi
            
            if total_ai < 0:
                total_ai = 0
            
            raw_ai_uptake += total_ai
            
            
            total_public_good = -1*secretor_public_good.uptakeInsideCellAtBoundaryTotalCount(cell, 0.5, 0.1).tot_amount
            
            
            # Make nutrient uptake modified by public good
            if total_public_good > 0:

                # Modify uptake values using modify function
                mod_nutrient_rel_uptake = modify_value(nutrient_rel_uptake, 2, total_public_good, exponential=False, linear=True)
                mod_nutrient_max_uptake = modify_value(nutrient_max_uptake, 2, total_public_good, exponential=False, linear=True)

                # Call total nutrient for both scenarios
                total_nutrient = -1*secretor_nutrient.uptakeOutsideCellAtBoundaryTotalCount(cell, 
                                                                                    mod_nutrient_max_uptake, 
                                                                                    mod_nutrient_rel_uptake).tot_amount
                
                sum_nutrient = total_nutrient
                
                if total_nutrient < 0.005:
                    # Create Chemotactic Preference for Complex
                    cd = self.chemotaxisPlugin.getChemotaxisData(cell, "complex_nutrient")
                    nd = self.chemotaxisPlugin.getChemotaxisData(cell, 'nutrient')
                    if cd:
                        # Increase complex nutrient chemotaxis
                        l = int(cd.getLambda() + 20.0)
                        cd.setLambda(l)
                        
                        # Decrease nutrient chemotaxis
                        n = int(nd.getLambda() / 200)
                        nd.setLambda(n)

                    # Define uptake of complex nutrients - set to 1/10th original nutrient uptake
                    complex_max_uptake = nutrient_max_uptake 
                    complex_rel_uptake = nutrient_rel_uptake 
                    
                    # Modify uptake values using modify function
                    mod_complex_rel = modify_value(complex_rel_uptake, 10, total_public_good)
                    mod_complex_max = modify_value(complex_rel_uptake, 10, total_public_good)
                    
                    # Call total complex nutrient
                    total_complex = -1 * secretor_complex.uptakeOutsideCellAtBoundaryTotalCount(cell, 
                                                                                                mod_complex_max, 
                                                                                                mod_complex_rel).tot_amount

                    sum_nutrient = total_nutrient + total_complex
                
            else:   
                total_nutrient = -1*secretor_nutrient.uptakeOutsideCellAtBoundaryTotalCount(cell, 
                                                                                            nutrient_max_uptake, 
                                                                                            nutrient_rel_uptake).tot_amount
                sum_nutrient = total_nutrient
            
            # Secretion of Autoinducer
            ai_secretion_factor =  cell.volume / cell.dict['Me']
            secretor_ai.secreteOutsideCellAtBoundary(cell, ai_secretion_factor)
            
            growth_factor = 0.001 * cell.dict['Me']
            
            # Growth Mutation of Cells - increases during starvation
            if not (mcs % 10):
                
                cell.targetVolume += growth_factor
                
                if mutation:
                    new_type_index = np.random.randint(1,4)
                    # if cell.dict['Me'] >= 10:
                        # if np.random.uniform(0.000, 1.000) < mutation_prob:
                            # cell.type = new_type_index
                    
                    if sum_nutrient < 1 and mcs > 500 and cell.dict['Me'] < 500:
                        if np.random.uniform(0.000, 1.000) < starvation_mutation_prob:
                            cell.type = new_type_index

             # WT Metabolic Modulations
            if cell.type == 1:
                wt_decay_modifier = 0
                if cell.dict['Me'] > 0:
                    
                    # Set ceiling of Me
                    # if cell.dict['Me'] > 2000:
                        # cell.dict['Me'] = 2000
                    
                    # Secretion of Public Good 
                    if not (mcs % 10):
                        # Update Me Plot
                        self.plot_me.add_data_point("WT Me", mcs, cell.dict['Me'])
                    # QS Loops
                    if cell.dict['sensor'] == 0 and total_ai > ai_threshold:
                        cell.dict['sensor'] = 1
                    
                    if cell.dict['sensor'] == 1 and total_ai < ai_lower_thres:
                        cell.dict['sensor'] = 0

                    # Secretion of Public Good
                    if cell.dict['sensor'] == 1:
                       if not qsi_active:
                        secretor_public_good.secreteOutsideCellAtBoundary(cell, 1)

                        
                       if qsi_active:
                           x = total_qsi
                           
                           if x > 0:
                                # IC50 of PTSP QSI is 0.35 uM - reciprocal function to describe this is y = 0.1 / (x - 0.15)
                                y = (0.1 / (x - 0.15))
                                
                                if y <=1:
                                    raw_public_good_excretion += secretor_public_good.secreteOutsideCellAtBoundary(cell, y)
                            
                           else:
                                raw_public_good_excretion += secretor_public_good.secreteOutsideCellAtBoundary(cell, 1)
                                
                       
                       if s_decay_modifier:
                           wt_decay_modifier += 0.0005
                       
                       # Create QS Starvation modifier
                       if starvation_bonus:
                           if cell.dict['Me'] <= 500 and sum_nutrient <= 0.05:
                               wt_decay_modifier /= 10
                               
                               # Modify growth
                               cell.targetVolume -= growth_factor
                               growth_factor /= 10
                               
                               cell.targetVolume += growth_factor
                               
                
                if s_decay_modifier and cell.dict['sensor'] == 1:
                    wt_decay_modifier += 0.0005
                
                cell.dict['Me'] = cell.dict['Me'] + total_nutrient - (decay_rate + wt_decay_modifier + growth_factor)
                    

            # UC Metabolic Modulations
            if cell.type == 2:
                uc_decay_modifier = 0
                cell.dict['sensor'] = 1
                if cell.dict['Me'] > 0:
                    
                    # Set ceiling of Me
                    # if cell.dict['Me'] > 2000:
                        # cell.dict['Me'] = 2000
                    
                    # Add decay penalty for QS
                    if s_decay_modifier and cell.dict['sensor'] == 1:
                        uc_decay_modifier += 0.0005
                        
                        if total_ai < ai_threshold:
                            dis_to_thres = ai_threshold - total_ai
                            
                            early_qs_penalty = 3 * np.log10(141.4 * dis_to_thres + 1)
                            
                            uc_decay_modifier += early_qs_penalty
                    
                    if not qsi_active:
                        secretor_public_good.secreteOutsideCellAtBoundary(cell, 1)
                        uc_decay_modifier += 1
                            
                    if qsi_active:
                        x = total_qsi
                               
                        if x > 0:
                            # IC50 of PTSP QSI is 0.35 uM - reciprocal function to describe this is y = 0.1 / (x - 0.15)
                            y = (0.1 / (x - 0.15))
                            if y <= 1:  
                                raw_public_good_excretion += secretor_public_good.secreteOutsideCellAtBoundary(cell, y)
                                uc_decay_modifier += 1
                        else:
                                    
                            raw_public_good_excretion += secretor_public_good.secreteOutsideCellAtBoundary(cell, 1)
                            uc_decay_modifier += 1
                    
                    # Secretion of Public Good and Mutation of WT - dice rolls every 10 steps
                    if not (mcs % 10):
                        # Update Me Plot
                        self.plot_me.add_data_point("UC Me", mcs, cell.dict['Me'])
                    
                cell.dict['Me'] = cell.dict['Me'] + total_nutrient - (decay_rate + uc_decay_modifier + growth_factor)
                
                
            
            # D Metabolic Modulations
            if cell.type == 3:
                d_decay_modifier = 0
                cell.dict['sensor'] = 0
                if cell.dict['Me'] > 0:
                    
                    # Set ceiling of Me
                    # if cell.dict['Me'] > 2000:
                        # cell.dict['Me'] = 2000
                    
                    # Set metabolic consequence of no qs genes - use ai field as measure of population divide by nutrients to ease if food is readily avaliable
                    if cheating_consequence:

                        # Set exponential decrease in metabolic consequence if nutrients are readily avaliable
                        
                        try:
                            metabolic_modifier = (total_ai / sum_nutrient)
                            if metabolic_modifier > 1:
                                metabolic_modifier = 1
                        
                        except:
                            metabolic_modifier = 100
                        
                        # Add to delta decay modifier
                        d_decay_modifier += metabolic_modifier
                    
                    cell.dict['Me'] = cell.dict['Me'] + total_nutrient - (decay_rate + d_decay_modifier + growth_factor)
                   
                if not (mcs % 10):
                    # Update Me Plot
                    
                    self.plot_me.add_data_point("D Me", mcs, cell.dict['Me'])  

        for cell in self.cell_list:
            # Plot New Data
            if not (mcs % 10):
                #self.plot_me.add_data_point("Me", mcs, cell.dict['Me'])
                
                wt_pop = len(self.cell_list_by_type(self.WT))
                uc_pop = len(self.cell_list_by_type(self.UC))
                d_pop = len(self.cell_list_by_type(self.D))
                
                tot_pop = wt_pop + uc_pop + d_pop
                
                self.plot_cell_pops.add_data_point('WT', mcs, (wt_pop/tot_pop) *  100)
                self.plot_cell_pops.add_data_point('UC', mcs, (uc_pop/tot_pop) * 100)
                self.plot_cell_pops.add_data_point('D', mcs, (d_pop/tot_pop) * 100)
                
                if track_cooperators:
                    if not (mcs % 100):
                        self.plot_coop.add_data_point('Cooperating Fraction', mcs, ((wt_pop + uc_pop) / tot_pop))
        
        # Calculate Avg AI Uptake
        avg_ai_uptake = raw_ai_uptake / len(self.cell_list)
        
        if qsi_active:
            if len(self.cell_list_by_type(self.WT)) + len(self.cell_list_by_type(self.UC)) > 0:
                avg_qsi_uptake = raw_qsi_uptake / (len(self.cell_list_by_type(self.WT)) + len(self.cell_list_by_type(self.UC)))
                
                avg_public_good_excretion = raw_public_good_excretion / (len(self.cell_list_by_type(self.WT)) + len(self.cell_list_by_type(self.UC)))
        
        if not (mcs % 10):
            self.plot_ai.add_data_point('AI Uptake', mcs, avg_ai_uptake)
            if qsi_active and len(self.cell_list_by_type(self.WT)) + len(self.cell_list_by_type(self.UC)) > 0:
                self.plot_qsi.add_data_point('QSI Uptake', mcs, avg_qsi_uptake)
                
                self.plot_qsi.add_data_point('Public Good Excretion', mcs, avg_public_good_excretion)
        # self.plot_cell_pops.add_data_point('WT', mcs, len(self.cell_list_by_type(self.WT))


 
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,frequency=1):
        MitosisSteppableBase.__init__(self,frequency)

    def step(self, mcs):

        cells_to_divide=[]
        for cell in self.cell_list:
            
            # Divisin gated by volume and metabolic energy
            if cell.volume > 100 and cell.dict['Me'] > 200:
                cells_to_divide.append(cell)

        for cell in cells_to_divide:

            self.divide_cell_random_orientation(cell)
            # Other valid options
            # self.divide_cell_orientation_vector_based(cell,1,1,0)
            # self.divide_cell_along_major_axis(cell)
            # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        # Division halves volume and metabolic energy
        if self.parent_cell.targetVolume > 100 and self.parent_cell.dict['Me'] > 200:
            
            self.parent_cell.targetVolume /= 2.0
            self.parent_cell.dict['Me'] /= 2.0
        

        self.clone_parent_2_child()            

   
class DeathSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        
        for cell in self.cell_list:
            # Kill cells randomly during dilution to mimic back dilution
            if dilution_bottleneck and dilution:
                if not (mcs % distance_btw_dilutions) and mcs > 0:
                    cell.dict['dilution_score'] = np.random.uniform(0.0, 1.0)
                
                    if cell.dict['dilution_score'] < dilution_factor:
                        cell.dict['Me'] = 0
            
            if cell.dict['Me'] <= 0:
                
                # Release Nutrients on Cell death
                # if (mcs % distance_btw_dilutions) > 0:
                    # secretor_nutrient = self.get_field_secretor("nutrient")
                    # secretor_nutrient.secreteOutsideCellAtBoundary(cell, 10)
                    
                    # secretor_complex = self.get_field_secretor('complex_nutrient')
                    # secretor_complex.secreteOutsideCellAtBoundary(cell, 20)
                    
                    # secretor_ai = self.get_field_secretor('ai')
                    # secretor_ai.secreteOutsideCellAtBoundary(cell, 10)
                
                # Cell Death
                cell.targetVolume=0
                cell.lambdaVolume=10000
                self.delete_cell(cell)


