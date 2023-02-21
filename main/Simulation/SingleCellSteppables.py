
from cc3d.core.PySteppables import *


class SingleCellSteppable(SteppableBasePy):

    def __init__(self, frequency=1):

        SteppableBasePy.__init__(self,frequency)

    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """
        x = self.dim.x//2
        y = self.dim.y//2
        size = 8
        cell = self.new_cell(self.CELL)
        # size of cell will be SIZExSIZEx1
        self.cell_field[x-size//2:x + size//2 - 1, 
                        y-size//2:y + size//2 - 1, 
                        0] = cell
        cell.lambdaVolume = 2.0
        cell.targetVolume = size * size
        
        
        
        
        
        

    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation
        
        :param mcs: current Monte Carlo step
        """

        for cell in self.cell_list:

            print("cell.id=",cell.id)

    def finish(self):
        """
        Called after the last MCS to wrap up the simulation
        """

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
