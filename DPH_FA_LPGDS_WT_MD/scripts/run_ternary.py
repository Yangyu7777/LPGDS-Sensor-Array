"""Two-phase safety/production runner for local binary-pair exploration."""
from pathlib import Path
import sys,json,time,os,csv
import numpy as np
from scipy.spatial import cKDTree
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
from openmm import app,unit,XmlSerializer,Platform,LangevinMiddleIntegrator,CustomExternalForce,MonteCarloBarostat
root=Path(__file__).resolve().parents[1];mode='FA';phase=sys.argv[1];replica=2
if phase not in ['smoke','production']:raise ValueError('Use smoke or production')
name='FA_R2';source=root/'inputs';out=root/'reproduced_ternary';statuspath=out/'status.json';seed=2026102102;ligand_name='DPH'
if phase=='smoke':
 if out.exists():raise RuntimeError('Output exists; refusing duplicate smoke run')
 out.mkdir();status={'status':'starting','mode':mode,'replica':replica,'independent_velocity_repeat':True,'source_endpoint':'Selected_MD_20261018/01_PROTEIN_DPH_R3, total5502ps','comparison':'Resident DPH; exterior competitor or matched control; no pushing force','pid':os.getpid(),'stage':'minimization','completed_stages':[],'planned_production_ps':5000,'planned_equilibration_ps':500,'short_test_ps':2,'timestep_fs':1,'random_seed':seed,'full_save_interval_ps':10,'solute_save_interval_ps':1,'limitations':'Independent velocity replica conditional on a fixed locally chosen start; no affinity or guaranteed binding; no experimental concentration match; provisional solvent model'}
else:
 status=json.loads(statuspath.read_text())
 if status['status']!='ready_for_production':raise RuntimeError('Safety stage not ready or production already launched')
 status.update(pid=os.getpid(),status='starting')
def save():
 temp=out/'status.tmp';temp.write_text(json.dumps(status,indent=2)+'\n');temp.replace(statuspath)
save()
try:
 p=app.PDBFile(str(source/'solvated.pdb'));p.positions=np.load(source/'coordinates_nm.npy')*unit.nanometer;atoms=list(p.topology.atoms());solute=[a.index for a in atoms if a.residue.name not in {'HOH','DMS','PH1','PH2','NA','K','CL'}];fa=[a.index for a in atoms if a.residue.name==ligand_name and a.element.symbol!='H'];partner=[a.index for a in atoms if a.index in solute and a.residue.name not in ['DPH','FAM','PFO'] and a.element.symbol!='H']
 if not fa or not partner:raise RuntimeError('Missing pair member')
 integrator=LangevinMiddleIntegrator(50*unit.kelvin,1/unit.picosecond,.001*unit.picosecond);integrator.setRandomNumberSeed(seed)
 if phase=='smoke':
  system=XmlSerializer.deserialize((source/'system.xml').read_text());restraint=CustomExternalForce('0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2');restraint.addGlobalParameter('k',1000)
  for n in ['x0','y0','z0']:restraint.addPerParticleParameter(n)
  px=np.asarray(p.positions.value_in_unit(unit.nanometer))
  for i in solute:
   if atoms[i].element.symbol!='H':restraint.addParticle(i,px[i].tolist())
  system.addForce(restraint);barostat=MonteCarloBarostat(1*unit.bar,298.15*unit.kelvin,0);barostat.setRandomNumberSeed(seed+51);system.addForce(barostat)
 else:
  system=XmlSerializer.deserialize((out/'running-system.xml').read_text());barostat=next(f for f in system.getForces() if isinstance(f,MonteCarloBarostat))
 backend=Platform.getPlatformByName('OpenCL');sim=app.Simulation(p.topology,system,integrator,backend,{'DeviceIndex':'0','Precision':'single'})
 if phase=='smoke':
  sim.context.setPositions(p.positions);sim.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=1000);sim.context.applyConstraints(1e-6);sim.context.setVelocitiesToTemperature(50*unit.kelvin,seed);sim.context.applyVelocityConstraints(1e-6)
 else:
  sim.loadCheckpoint(str(out/'latest.chk'));sim.currentStep=2000
  if abs(sim.context.getState().getTime().value_in_unit(unit.picosecond)-2)>1e-5:raise RuntimeError('Wrong safety checkpoint')
 status.update(status='running',device=backend.getPropertyValue(sim.context,'DeviceName'),precision='single');save()
 mols={};mapping={};torsions={}
 for n in {'CONTROL':['DPH'],'FA':['DPH','FAM'],'PFOA':['DPH','PFO']}[mode]:
  mol=Chem.SDMolSupplier(str(source/f'{n}.sdf'),removeHs=False)[0];mols[n]=mol;mapping[n]=[a.index for a in atoms if a.residue.name==n]
  if len(mapping[n])!=mol.GetNumAtoms():raise RuntimeError('Small molecule atom mismatch')
  ts=[]
  for b in mol.GetBonds():
   if b.GetBondType()!=Chem.BondType.DOUBLE or b.IsInRing() or b.GetBeginAtom().GetAtomicNum()!=6 or b.GetEndAtom().GetAtomicNum()!=6:continue
   i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx();aa=next(z.GetIdx() for z in mol.GetAtomWithIdx(i).GetNeighbors() if z.GetIdx()!=j and z.GetAtomicNum()!=1);bb=next(z.GetIdx() for z in mol.GetAtomWithIdx(j).GetNeighbors() if z.GetIdx()!=i and z.GetAtomicNum()!=1);ts.append((aa,i,j,bb))
  torsions[n]=ts
 def check():
  st=sim.context.getState(getPositions=True,getEnergy=True);x=np.asarray(st.getPositions().value_in_unit(unit.angstrom));energy=st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
  if not np.isfinite(x).all() or not np.isfinite(energy):raise RuntimeError('Non-finite state')
  angles={}
  for n,mol in mols.items():
   for i,j in enumerate(mapping[n]):mol.GetConformer().SetAtomPosition(i,x[j])
   angles[n]=[rdMolTransforms.GetDihedralDeg(mol.GetConformer(),*t) for t in torsions[n]]
   if any(abs(z)<90 for z in angles[n]):raise RuntimeError(f'{n} E/Z geometry changed: {angles[n]}')
  box=np.diag(np.asarray(st.getPeriodicBoxVectors().value_in_unit(unit.angstrom)));tree=cKDTree(x[partner]%box,boxsize=box);gap=float(tree.query(x[fa]%box)[0].min())
  status.update(step=sim.currentStep,time_ps=float(st.getTime().value_in_unit(unit.picosecond)),min_ligand_partner_heavy_A=gap,small_molecule_dihedrals_deg=angles,potential_energy_kj_mol=energy,restraint_k=sim.context.getParameter('k'));sim.saveCheckpoint(str(out/'latest.chk'));sim.saveState(str(out/'latest-state.xml'));save()
 if phase=='smoke':
  status['stage']='short_test';save();sim.step(2000);check()
  (out/'running-system.xml').write_text(XmlSerializer.serialize(system));status.update(status='ready_for_production',short_test_passed=True,minimization_passed=True);save();print(mode,'passed minimization and 2 ps safety test',flush=True)
 else:
  sim.reporters.append(app.StateDataReporter(str(out/'thermodynamics.csv'),1000,step=True,time=True,potentialEnergy=True,kineticEnergy=True,temperature=True,volume=True,density=True,speed=True));sim.reporters.append(app.DCDReporter(str(out/'solute-trajectory.dcd'),1000,atomSubset=solute,enforcePeriodicBox=False));sim.reporters.append(app.DCDReporter(str(out/'full-trajectory.dcd'),10000,enforcePeriodicBox=False))
  with (out/'full-topology.pdb').open('w') as f:app.PDBFile.writeFile(p.topology,p.positions,f)
  sm=app.Modeller(p.topology,p.positions);ss=set(solute);sm.delete([a for a in atoms if a.index not in ss])
  with (out/'solute-topology.pdb').open('w') as f:app.PDBFile.writeFile(sm.topology,sm.positions,f)
  stages=[('heating',20,1000,False),('NVT',80,1000,False),('NPT_strong',100,1000,True),('NPT_weak',100,100,True),('NPT_release',200,0,True),('production_pilot',5000,0,True)];start=time.monotonic()
  for stage,ps,k,pressure in stages:
   status.update(stage=stage,stage_duration_ps=ps);save();sim.context.setParameter('k',k);barostat.setFrequency(25 if pressure else 0);sim.context.reinitialize(preserveState=True);print(mode,'starting',stage,flush=True)
   for block in range(ps//5):
    temperature=50+(298.15-50)*(block+1)/(ps//5) if stage=='heating' else 298.15;integrator.setTemperature(temperature*unit.kelvin);sim.step(5000);check()
   if stage=='NPT_release':
    with (out/'thermodynamics.csv').open() as f:data=list(csv.reader(f))
    header=[z.strip('#\"') for z in data[0]];recent=data[-100:];temps=np.array([float(z[header.index('Temperature (K)')]) for z in recent]);density=np.array([float(z[header.index('Density (g/mL)')]) for z in recent]);status['preproduction_screen']={'mean_temperature_K':float(temps.mean()),'mean_density_g_ml':float(density.mean()),'density_range_g_ml':float(np.ptp(density))}
    if not (285<temps.mean()<310 and .85<density.mean()<1.2 and np.ptp(density)<.1):raise RuntimeError('Preproduction sanity screen failed')
   st=sim.context.getState(getPositions=True);p.topology.setPeriodicBoxVectors(st.getPeriodicBoxVectors())
   with (out/f'{stage}-final.pdb').open('w') as f:app.PDBFile.writeFile(p.topology,st.getPositions(),f)
   (out/'running-system.xml').write_text(XmlSerializer.serialize(system));status['completed_stages'].append(stage);save()
  status.update(status='complete',elapsed_seconds=time.monotonic()-start,production_completed_ps=5000);save()
except BaseException as e:status.update(status='failed',error=str(e));save();raise
