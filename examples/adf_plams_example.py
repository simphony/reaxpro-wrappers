from scm.plams import *
init()
sett = Settings()
sett.input.ams.task='GeometryOptimization'
sett.input.ADF.basis.type = 'DZP'
sett.input.ADF.xc.gga = 'PBE'

mol = Molecule()
mol.add_atom(Atom(symbol='O', coords=(0,0,0)))
mol.add_atom(Atom(symbol='H', coords=(1,0,0)))
mol.add_atom(Atom(symbol='H', coords=(0,1,0)))

# Alternatively you can define the molecule loading the geomtry from
# a .xyz file:
# mol = Molecule('XYZ/H2O.xyz')

# Run ADF job:
job = AMSJob(molecule=mol, settings=sett)

# Run job:
results = job.run()
# Print results:
nEntries = results.readrkf('History', 'nEntries')
print("Total energy per geometry optimization step:")
for k in range(nEntries):
    s = "Energy({})".format(k+1)
    x = results.readrkf("History", s)
    print("Optimization cycle", k+1, ":", x, "Hartree")
