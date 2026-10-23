from pathlib import Path
import sys,json,csv
import numpy as np
from rdkit import Chem
from openmm import app
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'));from md_analysis_utils import read_dcd
x,b=read_dcd(R/'raw/solute-trajectory.dcd');at=list(app.PDBFile(str(R/'raw/solute-topology.pdb')).topology.atoms());t=np.arange(3,5503);data=[];names=[];stats={}
def dihedral(q):
 a=-(q[:,1]-q[:,0]);v=q[:,2]-q[:,1];v/=np.linalg.norm(v,axis=1)[:,None];c=q[:,3]-q[:,2];aa=a-(a*v).sum(axis=1)[:,None]*v;cc=c-(c*v).sum(axis=1)[:,None]*v;return np.degrees(np.arctan2((np.cross(v,aa)*cc).sum(axis=1),(aa*cc).sum(axis=1)))
for n,expect in [('DPH','E'),('FAM','Z')]:
 rd=Chem.SDMolSupplier(str(R/'inputs'/f'{n}.sdf'),removeHs=False)[0];idx=[a.index for a in at if a.residue.name==n];angles=[]
 for bd in rd.GetBonds():
  if bd.GetBondType()!=Chem.BondType.DOUBLE or bd.IsInRing() or bd.GetBeginAtom().GetAtomicNum()!=6 or bd.GetEndAtom().GetAtomicNum()!=6:continue
  i,j=bd.GetBeginAtomIdx(),bd.GetEndAtomIdx();a=next(z.GetIdx() for z in rd.GetAtomWithIdx(i).GetNeighbors() if z.GetIdx()!=j and z.GetAtomicNum()!=1);c=next(z.GetIdx() for z in rd.GetAtomWithIdx(j).GetNeighbors() if z.GetIdx()!=i and z.GetAtomicNum()!=1);ang=dihedral(x[:,[idx[k] for k in [a,i,j,c]]].astype(float));angles.append(ang);data.append(ang);names.append(f'{n}_CC_double_bond_{len(angles)}_dihedral_deg')
 aa=np.array(angles);good=(np.abs(aa)>90).all() if expect=='E' else (np.abs(aa)<90).all();stats[n]={'expected':expect,'double_bonds':len(angles),'all_5500_saved_frames_correct':bool(good),'min_absolute_dihedral_deg':float(np.abs(aa).min()),'max_absolute_dihedral_deg':float(np.abs(aa).max())};assert good
np.savetxt(R/'analysis/double_bond_dihedrals.csv',np.c_[t,*data],delimiter=',',header=','.join(['total_time_ps']+names),comments='');(R/'provenance/stereochemistry_all_saved_frames.json').write_text(json.dumps(stats,indent=2)+'\n');print(json.dumps(stats,indent=2))
