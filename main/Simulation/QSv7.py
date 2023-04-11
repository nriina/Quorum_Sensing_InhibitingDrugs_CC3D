
from cc3d import CompuCellSetup
        


from QSv7Steppables import ConstraintInitializerSteppable

CompuCellSetup.register_steppable(steppable=ConstraintInitializerSteppable(frequency=1))




from QSv7Steppables import GrowthSteppable

CompuCellSetup.register_steppable(steppable=GrowthSteppable(frequency=1))




from QSv7Steppables import MitosisSteppable

CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))




from QSv7Steppables import DeathSteppable

CompuCellSetup.register_steppable(steppable=DeathSteppable(frequency=1))


CompuCellSetup.run()
