from pathlib import Path
import sys,json,csv
import numpy as np
from openmm import app,unit
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'));from md_analysis_utils import read_dcd,kabsch
m=np.genfromtxt(R/'analysis/metrics.csv',delimiter=',',names=True);joint=(m['DPH_protein_min_heavy_A']<4.5)&(m['competitor_protein_min_heavy_A']<4.5)&(m['DPH_competitor_min_heavy_A']<4.5)
edges=np.diff(np.r_[False,joint,False].astype(int));starts=np.where(edges==1)[0];ends=np.where(edges==-1)[0];imax=np.argmax(ends-starts);k=(starts[imax]+ends[imax]-1)//2;t=int(m['total_time_ps'][k]);top=app.PDBFile(str(R/'raw/solute-topology.pdb'));at=list(top.topology.atoms());ref=np.asarray(app.PDBFile(str(R/'snapshots/aligned_snapshot_0000ps.pdb')).positions.value_in_unit(unit.angstrom));x,b=read_dcd(R/'raw/solute-trajectory.dcd');q=x[k+500].astype(float);box=b[k+500];fit=np.array([a.index for a in at if a.name=='CA' and a.residue.name not in ['DPH','FAM']])[:170];rot,mc,tc=kabsch(q[fit],ref[fit]);aligned=(q-mc)@rot+tc
for n in ['DPH','FAM']:
 sel=[a.index for a in at if a.residue.name==n];shift=np.round((q[sel].mean(axis=0)-mc)/box)*box;aligned[sel]=(q[sel]-shift-mc)@rot+tc
top.topology.setPeriodicBoxVectors(None);fname=f'aligned_representative_contact_{t:04d}ps.pdb'
with (R/'snapshots'/fname).open('w') as f:app.PDBFile.writeFile(top.topology,aligned*unit.angstrom,f)
j={'selection_rule':'Midpoint saved frame of the longest contiguous episode where all three pairwise minimum heavy-atom distances are <4.5 A. This is a selected contact illustration, not an independently representative equilibrium state.','selected_total_time_ps':t,'selected_production_time_ns':float(m['production_time_ns'][k]),'episode_start_production_ns':float(m['production_time_ns'][starts[imax]]),'episode_last_production_ns':float(m['production_time_ns'][ends[imax]-1]),'episode_saved_frames':int(ends[imax]-starts[imax]),'episode_estimated_duration_ns':int(ends[imax]-starts[imax])/1000,'snapshot':fname,'min_distances_A':{n:float(m[n][k]) for n in ['DPH_protein_min_heavy_A','competitor_protein_min_heavy_A','DPH_competitor_min_heavy_A']}}
(R/'provenance/contact_snapshot_selection.json').write_text(json.dumps(j,indent=2)+'\n')
with (R/'analysis/three_pair_contacts.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['production_time_ns','DPH_WT_contact','C18_2_WT_contact','DPH_C18_2_contact','all_three_pair_contacts']);w.writerows(zip(m['production_time_ns'],m['DPH_protein_min_heavy_A']<4.5,m['competitor_protein_min_heavy_A']<4.5,m['DPH_competitor_min_heavy_A']<4.5,joint))
with (R/'analysis/contact_threshold_sensitivity.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['cutoff_A','DPH_WT_fraction','FA_WT_fraction','DPH_FA_fraction','all_three_fraction'])
 for cut in [3.5,4.,4.5,5.]:
  a=m['DPH_protein_min_heavy_A']<cut;b=m['competitor_protein_min_heavy_A']<cut;c=m['DPH_competitor_min_heavy_A']<cut;w.writerow([cut,a.mean(),b.mean(),c.mean(),(a&b&c).mean()])
print(json.dumps(j,indent=2))
