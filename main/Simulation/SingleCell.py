
from cc3d import CompuCellSetup
        

from SingleCellSteppables import SingleCellSteppable

CompuCellSetup.register_steppable(steppable=SingleCellSteppable(frequency=1))


CompuCellSetup.run()
