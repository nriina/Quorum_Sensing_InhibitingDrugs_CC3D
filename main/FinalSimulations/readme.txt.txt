2/19 Logan and I started a base simulation

cells have 3 dictionary attributes

sensor = sensor allele
sensornum = sensor id for tracking field
ligand = ligand allele
Me = metabolic energy

we have implemented:

metabolic energy that decays according to # of neighbors that have S allele
growth
death if me <= 0

tracking fields for Me, and sensornum (sensor allele)



# to do next
make a discintion between cells locked in S state and cells that switch to S-state based on chemical field
once we go over chemical fields in class, implement chemical fields
think about how ligand allele will impact Me
	producing will have a metabolic cost,
	there will be a production rate which will be dependent on chemical field, which will effect metabolic energy
