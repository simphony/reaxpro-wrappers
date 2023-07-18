from scm.plams import *
init()
sett = Settings()
sett.input.ams.task='GeometryOptimization'
sett.input.ADF.basis.type = 'TZP'
sett.input.ADF.xc.gga = 'PBE'

sett2= Settings()
sett2.input.ams.task='GeometryOptimization'
sett2.input.ADF.basis.type = 'DZP'
sett2.input.ADF.xc.gga = 'PBE'
mol2 = Molecule( './XYZ/H2O.xyz' )

mol = Molecule()
mol.add_atom(Atom(symbol='O', coords=(0,0,0)))
mol.add_atom(Atom(symbol='H', coords=(1,0,0)))
mol.add_atom(Atom(symbol='H', coords=(0,1,0)))

# Alternatively you can define the molecule loading the geomtry from
# a .xyz file:
# mol = Molecule('XYZ/H2O.xyz')

# Run ADF job:
job1 = AMSJob(molecule=mol, settings=sett)
job2 = AMSJob(molecule=mol, settings=sett2)
job_list=[job1, job2]
multijob=MultiJob()
multijob.children=job_list

# Run job:
#results = job.run()
#results_multijob = multijob.run()
multijob.run()
job2E = multijob.children[1].results.get_energy(unit='eV')
print("Total energy:", job2E, "eV") 
## Print results:
#nEntries = results.readrkf('History', 'nEntries')
#print("Total energy per geometry optimization step:")
#for k in range(nEntries):
#    s = "Energy({})".format(k+1)
#    x = results.readrkf("History", s)
#    print("Optimization cycle", k+1, ":", x, "Hartree")
