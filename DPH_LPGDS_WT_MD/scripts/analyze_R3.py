"""Reanalyse preserved trajectories; no simulation is launched. Units: angstrom, ps."""
from pathlib import Path
import sys,json,csv,hashlib,platform,collections
import numpy as np
from scipy.spatial import cKDTree
from openmm import app,unit,XmlSerializer,NonbondedForce
from md_analysis_utils import read_dcd,kabsch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.titlesize':11})
COLORS={'DPH':'#16657a','FAM':'#a75a32'}
def writecsv(path,header,rows):
 with path.open('w') as f:w=csv.writer(f);w.writerow(header);w.writerows(rows)
def analyze(folder,ligand):
 d=BASE;out=d/'analysis';figdir=d/'figures';raw=d/'raw';src=d/'inputs';col=COLORS[ligand];label='Protein–DPH R3' if ligand=='DPH' else 'Protein–FA (20261015)'
 p=app.PDBFile(str(raw/'solute-topology.pdb'));at=list(p.topology.atoms());xyz,boxes=read_dcd(raw/'solute-trajectory.dcd');time=np.arange(3,len(xyz)+3,dtype=float);mask=time>502;t=(time[mask]-502)/1000;assert mask.sum()==5000 and time[-1]==5502
 lig=np.array([a.index for a in at if a.residue.name==ligand]);lh=np.array([i for i in lig if at[i].element.symbol!='H']);prot=np.array([a.index for a in at if a.residue.name!=ligand]);ph=np.array([i for i in prot if at[i].element.symbol!='H']);ca=np.array([i for i in prot if at[i].name=='CA']);fit=ca[:170];res=list(p.topology.residues())[:-1];resmap={r.index:j for j,r in enumerate(res)};labels=[r.name+r.id for r in res]
 initial=app.PDBFile(str(src/'solvated.pdb'));ia=list(initial.topology.atoms());ix=np.asarray(initial.positions.value_in_unit(unit.angstrom));keep=[a.index for a in ia if a.residue.name not in {'HOH','DMS','PH1','PH2','NA','K','CL'}];ref=ix[keep];assert len(ref)==len(at)
 mass=np.array([at[i].element.mass.value_in_unit(unit.dalton) for i in ph]);refcen=ref[lh].mean(axis=0)
 contact=np.zeros((len(xyz),len(res)),bool);metrics=[];ca_all=[];lig_all=[];snapshots={0:ref.copy()};wanted={502,1502,3002,5502}
 for k,(x,box,tm) in enumerate(zip(xyz,boxes,time)):
  x=x.astype(float);rot,mc,tc=kabsch(x[fit],ref[fit]);aligned=(x-mc)@rot+tc;shift=np.round((x[lh].mean(axis=0)-mc)/box)*box;aligned[lig]=(x[lig]-shift-mc)@rot+tc
  tree=cKDTree(x[ph]%box,boxsize=box);dist=tree.query(x[lh]%box)[0];hits=tree.query_ball_point(x[lh]%box,4.5)
  for j in set(z for h in hits for z in h):contact[k,resmap[at[ph[j]].residue.index]]=True
  prms=np.sqrt(np.mean(np.sum((aligned[fit]-ref[fit])**2,axis=1)));pcen=np.average(aligned[ph],axis=0,weights=mass);rg=np.sqrt(np.average(np.sum((aligned[ph]-pcen)**2,axis=1),weights=mass));lfit=np.sqrt(np.mean(np.sum((aligned[lh]-ref[lh])**2,axis=1)));r,m,c=kabsch(aligned[lh],ref[lh]);selfrms=np.sqrt(np.mean(np.sum(((aligned[lh]-m)@r+c-ref[lh])**2,axis=1)));lc=aligned[lh].mean(axis=0)
  metrics.append([tm,(tm-502)/1000,dist.min(),contact[k].sum(),prms,rg,lfit,selfrms,np.linalg.norm(lc-refcen),*lc]);ca_all.append(aligned[ca]);lig_all.append(aligned[lh])
  if tm in wanted:snapshots[int(tm)]=aligned.copy()
 metrics=np.asarray(metrics);ca_all=np.asarray(ca_all);lig_all=np.asarray(lig_all)
 header=['total_time_ps','production_time_ns','min_heavy_distance_A','contact_residue_count','protein_core_CA_RMSD_preMD_A','protein_heavy_mass_weighted_Rg_A','ligand_protein_aligned_RMSD_preMD_A','ligand_self_aligned_RMSD_preMD_A','ligand_centroid_displacement_preMD_A','ligand_centroid_x_A','ligand_centroid_y_A','ligand_centroid_z_A'];writecsv(out/'trajectory_metrics.csv',header,metrics[mask]);np.savez_compressed(out/'residue_contacts.npz',total_time_ps=time[mask],residue_labels=labels,contact=contact[mask]);occ=contact[mask].mean(axis=0);writecsv(out/'residue_contact_occupancy.csv',['residue','contact_frames','total_frames','occupancy_fraction'],[[l,int(contact[mask,j].sum()),5000,o] for j,(l,o) in enumerate(zip(labels,occ))])
 rf=np.sqrt(np.mean(np.sum((ca_all[mask]-ca_all[mask].mean(axis=0))**2,axis=2),axis=0));writecsv(out/'protein_CA_RMSF.csv',['residue','input_residue_id','RMSF_A'],[[at[i].residue.name+at[i].residue.id,at[i].residue.id,v] for i,v in zip(ca,rf)])
 threshold=[(v,int((metrics[mask,2]<v).sum()),float((metrics[mask,2]<v).mean())) for v in [3.5,4.0,4.5,5.0]];writecsv(out/'contact_threshold_sensitivity.csv',['cutoff_A','contact_frames','fraction_of_5000_frames'],threshold)
 blocks=[]
 for b in range(5):
  m=(t>b)&(t<=b+1);z=metrics[mask][m];blocks.append([b+1,len(z),np.mean(z[:,2]<4.5),*z[:,2:9].mean(axis=0)])
 writecsv(out/'one_ns_block_descriptive_means.csv',['block','frames','contact_fraction']+header[2:9],blocks)
 # Geometric ligand SASA with and without protein, excluding water/cosolvent as occluders.
 rad=np.array([{'H':1.2,'C':1.7,'N':1.55,'O':1.52,'S':1.8}[a.element.symbol] for a in at]);probe=1.4
 def sasa(frame,box,n=960):
  ii=np.arange(n);zz=1-2*(ii+.5)/n;phi=ii*np.pi*(3-np.sqrt(5));rr=np.sqrt(1-zz*zz);dirs=np.column_stack([rr*np.cos(phi),rr*np.sin(phi),zz]);lr=rad[lig]+probe;pr=rad[prot]+probe;iso=comp=0.
  for j,atom in enumerate(lig):
   pts=frame[atom]+dirs*lr[j];delta=pts[:,None,:]-frame[lig][None,:,:];delta-=np.round(delta/box)*box;blocked=np.sum(delta*delta,axis=2)<lr[None,:]**2-1e-8;blocked[:,j]=False;free=~blocked.any(axis=1);delta=frame[prot]-frame[atom];delta-=np.round(delta/box)*box;near=np.flatnonzero(np.linalg.norm(delta,axis=1)<lr[j]+pr);q=pts[:,None,:]-(frame[atom]+delta[near])[None,:,:];occluded=(np.sum(q*q,axis=2)<pr[near][None,:]**2-1e-8).any(axis=1);weight=4*np.pi*lr[j]**2/n;iso+=free.sum()*weight;comp+=(free&~occluded).sum()*weight
  return [iso,comp,1-comp/iso]
 ss=np.array([[time[k],(time[k]-502)/1000,*sasa(xyz[k].astype(float),boxes[k])] for k in np.flatnonzero(mask&(time%50==0))]);final=sasa(xyz[-1].astype(float),boxes[-1]);fine=sasa(xyz[-1].astype(float),boxes[-1],3840);writecsv(out/'ligand_SASA.csv',['total_time_ps','production_time_ns','isolated_SASA_A2','complex_SASA_A2','protein_occluded_fraction'],ss)
 print(label,'metrics and SASA calculated',flush=True)
 if ligand=='DPH':sv=np.genfromtxt(d/'legacy_analysis/solvent.csv',delimiter=',',names=True);water=sv['ligand_water_O_within5A'];dmso=sv['ligand_DMSO_within4_5A']
 else:sv=np.genfromtxt(d/'legacy_analysis/solvent-neighbors.csv',delimiter=',',names=True);water=sv['FA_water_O_within5A'];dmso=sv['FA_DMSO_molecules_within4_5A']
 writecsv(out/'solvent_neighbors.csv',['total_time_ps','production_time_ns','water_oxygen_within5A','DMSO_molecules_within4_5A'],zip(sv['time_ps'],(sv['time_ps']-502)/1000,water,dmso))
 td=np.genfromtxt(raw/'thermodynamics.csv',delimiter=',',names=True);thermonames=td.dtype.names
 # NumPy strips spaces and punctuation from StateDataReporter column names.
 temp=td[next(n for n in thermonames if n.startswith('Temperature'))];density=td[next(n for n in thermonames if n.startswith('Density'))];vol=td[next(n for n in thermonames if n.startswith('Box_Volume'))]
 assert len(td)==len(time)
 writecsv(out/'thermodynamics_production.csv',['total_time_ps','production_time_ns','temperature_K','density_g_mL','box_volume_nm3'],zip(time[mask],t,temp[mask],density[mask],vol[mask]))
 # Verify the endpoint independently from the full serialized state (not DCD).
 state=XmlSerializer.deserialize((raw/'latest-state.xml').read_text());sx=np.asarray(state.getPositions().value_in_unit(unit.angstrom));sb=np.diag(np.asarray(state.getPeriodicBoxVectors().value_in_unit(unit.angstrom)));sl=[a.index for a in ia if a.residue.name==ligand and a.element.symbol!='H'];sp=[keep[i] for i in ph];endpoint=cKDTree(sx[sp]%sb,boxsize=sb).query(sx[sl]%sb)[0].min();assert abs(endpoint-metrics[-1,2])<1e-4
 # Check initial packing distance, actual molecule counts, charge and PME settings.
 system=XmlSerializer.deserialize((src/'system.xml').read_text());nb=next(f for f in system.getForces() if isinstance(f,NonbondedForce));counts=collections.Counter(r.name for r in initial.topology.residues());initial_box=np.diag(np.asarray(initial.topology.getPeriodicBoxVectors().value_in_unit(unit.angstrom)));initial_volume=np.prod(initial_box)/1000
 initial_dist=float(cKDTree(ref[ph]%initial_box,boxsize=initial_box).query(ref[lh]%initial_box)[0].min())
 settings={'total_atoms':len(ia),'protein_residues':len(res),'core_residues':170,'His_tag_residues':6,'counts_solvent_ions':{n:counts[n] for n in ['HOH','DMS','PH1','PH2','NA','K','CL']},'initial_box_A':initial_box.tolist(),'initial_volume_nm3':initial_volume,'initial_protein_and_ligand_concentration_mM':1000/(initial_volume*.602214076),'production_mean_inverse_volume_protein_mM':float(np.mean(1000/(vol[mask]*.602214076))),'net_charge_e':float(sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(system.getNumParticles()))),'constraints':system.getNumConstraints(),'nonbonded_method_enum':int(nb.getNonbondedMethod()),'cutoff_nm':nb.getCutoffDistance().value_in_unit(unit.nanometer),'ewald_tolerance':nb.getEwaldErrorTolerance(),'switching':nb.getUseSwitchingFunction(),'dispersion_correction':nb.getUseDispersionCorrection(),'initial_min_heavy_distance_A':initial_dist,'hydrogen_masses_Da':sorted(set(round(system.getParticleMass(a.index).value_in_unit(unit.dalton),6) for a in ia if a.element.symbol=='H'))}
 summary={'system':label,'production_frames':5000,'production_duration_ns':5,'analysis_reference':'Pre-MD prepared structure; fit 170 core C-alpha atoms; His6 excluded from fit','contact_fraction':float((metrics[mask,2]<4.5).mean()),'contact_frames':int((metrics[mask,2]<4.5).sum()),'min_distance_A':float(metrics[mask,2].min()),'final_min_distance_A':float(metrics[-1,2]),'mean_core_CA_RMSD_A':float(metrics[mask,4].mean()),'mean_heavy_mass_weighted_Rg_A':float(metrics[mask,5].mean()),'mean_core_CA_RMSF_A':float(rf[:170].mean()),'ligand_final_protein_aligned_RMSD_A':float(metrics[-1,6]),'ligand_final_self_aligned_RMSD_A':float(metrics[-1,7]),'ligand_final_net_displacement_A':float(metrics[-1,8]),'mean_SASA_A2':float(ss[:,3].mean()),'final_SASA_A2':float(final[1]),'mean_protein_occluded_fraction':float(ss[:,4].mean()),'final_protein_occluded_fraction':float(final[2]),'final_SASA_resolution_difference_percentage_points':float(abs(final[2]-fine[2])*100),'SASA_frames':len(ss),'mean_water_O_within5A':float(water.mean()),'mean_DMSO_within4_5A':float(dmso.mean()),'mean_temperature_K':float(temp[mask].mean()),'mean_density_g_mL':float(density[mask].mean()),'top_contacts':[(labels[j],float(occ[j])) for j in np.argsort(-occ)[:15]],'independent_state_endpoint_A':float(endpoint),'endpoint_DCD_difference_A':float(abs(endpoint-metrics[-1,2])),'actual_system_settings':settings,'statistical_scope':'One deposited 5 ns trajectory; temporally correlated frame statistics are descriptive and do not establish binding affinity, spontaneous association, or enclosed-pocket residence.'}
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 def save(fig,name):
  for ext in ['png','pdf','svg']:fig.savefig(figdir/(name+'.'+ext),dpi=350,bbox_inches='tight')
  plt.close(fig)
 def style(ax,y):ax.set(xlabel='Production time (ns)',ylabel=y,xlim=(0,5));ax.grid(False)
 fig,axs=plt.subplots(3,1,figsize=(7.2,7.8),layout='constrained')
 for ax,y,v in [(axs[0],'Minimum heavy-atom distance (Å)',metrics[mask,2]),(axs[1],'Contacting residues',metrics[mask,3]),(axs[2],'Protein core Cα RMSD (Å)',metrics[mask,4])]:ax.plot(t,v,color=col,lw=.7);style(ax,y)
 axs[0].axhline(4.5,color='#555555',ls='--',lw=1,label='4.5 Å contact cutoff');axs[0].legend(frameon=False);fig.suptitle(label+' | Proximity and protein motion');save(fig,'Fig01_proximity_and_RMSD')
 fig,axs=plt.subplots(3,1,figsize=(7.2,7.6),layout='constrained');axs[0].plot(t,metrics[mask,5],color=col,lw=.7);style(axs[0],'Protein heavy-atom Rg (Å)');axs[1].plot(range(1,171),rf[:170],color=col);axs[1].set(xlabel='Input residue number',ylabel='Core Cα RMSF (Å)',xlim=(1,170));axs[2].plot(range(171,177),rf[170:],color=col,marker='o');axs[2].set(xlabel='Input residue number (His6 tag)',ylabel='Tag Cα RMSF (Å)',xlim=(171,176));fig.suptitle(label+' | Compactness and residue fluctuations');save(fig,'Fig02_Rg_and_RMSF')
 rank=np.argsort(-occ)[:12];matrix=contact[mask].reshape(50,100,len(res)).mean(axis=1)[:,rank].T*100;fig,(ax,b)=plt.subplots(1,2,figsize=(10,5.5),layout='constrained',gridspec_kw={'width_ratios':[1.5,1]});im=ax.imshow(matrix,aspect='auto',extent=[0,5,11.5,-.5],vmin=0,vmax=100,cmap='Blues');ax.set_yticks(range(12),[labels[j] for j in rank]);ax.set(xlabel='Production time (ns)',title='100 ps blocks');fig.colorbar(im,ax=ax,label='Contact frames (%)',shrink=.7);b.barh(range(12),occ[rank]*100,color=col);b.set_yticks(range(12),[labels[j] for j in rank]);b.invert_yaxis();b.set(xlabel='Contact occupancy (%)',xlim=(0,113));
 for i,j in enumerate(rank):b.text(occ[j]*100+1,i,f'{100*occ[j]:.1f}',va='center',fontsize=9)
 fig.suptitle(label+' | Residue contacts (<4.5 Å)');save(fig,'Fig03_residue_contacts')
 fig,axs=plt.subplots(3,1,figsize=(7.2,7.8),layout='constrained');axs[0].plot(ss[:,1],ss[:,2],color='#777777',ls='--',label='Isolated ligand, same conformation');axs[0].plot(ss[:,1],ss[:,3],color=col,label='Ligand + protein');style(axs[0],'Ligand SASA (Å²)');axs[0].legend(frameon=False,fontsize=9);axs[1].plot(ss[:,1],ss[:,4]*100,color=col);style(axs[1],'Protein-occluded SASA (%)');axs[1].set_ylim(0,100);axs[2].plot((sv['time_ps']-502)/1000,water,color=col);style(axs[2],'Water O atoms within 5 Å');fig.suptitle(label+' | Surface exposure and hydration');save(fig,'Fig04_SASA_and_hydration')
 fig,axs=plt.subplots(3,1,figsize=(7.2,7.8),layout='constrained');axs[0].plot(t,metrics[mask,6],color=col,label='Protein aligned');axs[0].plot(t,metrics[mask,7],color='#777777',ls='--',label='Ligand self aligned');style(axs[0],'Ligand RMSD (Å)');axs[0].legend(frameon=False);axs[1].plot(t,metrics[mask,8],color=col);style(axs[1],'Centroid displacement (Å)');axs[2].step((sv['time_ps']-502)/1000,dmso,color=col,where='mid');style(axs[2],'DMSO molecules within 4.5 Å');fig.suptitle(label+' | Ligand motion and cosolvent');save(fig,'Fig05_ligand_motion_and_DMSO')
 fig,axs=plt.subplots(2,1,figsize=(7.2,5.7),layout='constrained');axs[0].plot(t,temp[mask],color=col,lw=.6);style(axs[0],'Temperature (K)');axs[0].axhline(298.15,color='#777777',ls='--');axs[1].plot(t,density[mask],color=col,lw=.6);style(axs[1],'Density (g mL⁻¹)');fig.suptitle(label+' | Thermodynamic diagnostics');save(fig,'FigS01_temperature_and_density')
 # Common frame, common scale; actual snapshots, no artificial interpolation.
 p.topology.setPeriodicBoxVectors(None);bonds=[(a.index,b.index) for a,b in p.topology.bonds() if a.index in lh and b.index in lh];fig=plt.figure(figsize=(10.5,8),layout='constrained');frames=[snapshots[k] for k in [0,502,3002,5502]];points=np.concatenate([f[np.r_[ca,lh]] for f in frames]);center=(points.min(axis=0)+points.max(axis=0))/2;span=np.ptp(points,axis=0).max()/2+2
 for n,(tm,f) in enumerate(zip([0,502,3002,5502],frames),1):
  with (d/'snapshots'/f'aligned_{tm:04d}ps.pdb').open('w') as handle:app.PDBFile.writeFile(p.topology,f*unit.angstrom,handle)
  ax=fig.add_subplot(2,2,n,projection='3d');ax.plot(f[ca[:170],0],f[ca[:170],1],f[ca[:170],2],color='#a0a9b0',lw=1)
  for u,v in bonds:ax.plot(f[[u,v],0],f[[u,v],1],f[[u,v],2],color=col,lw=2.5)
  ax.set(xlim=(center[0]-span,center[0]+span),ylim=(center[1]-span,center[1]+span),zlim=(center[2]-span,center[2]+span),title='Pre-MD prepared structure' if tm==0 else f'Production {(tm-502)/1000:.1f} ns');ax.set_box_aspect((1,1,1));ax.view_init(elev=18,azim=-65);ax.set_axis_off()
 fig.suptitle(label+' | Aligned coordinate snapshots');save(fig,'FigS02_coordinate_snapshots')
 # PyMOL source for molecular surface figures, based on actual aligned coordinates.
 (d/'snapshots/view_in_pymol.pml').write_text('reinitialize\nload aligned_5502ps.pdb, final\nhide everything\nshow cartoon, polymer.protein\ncolor gray70, polymer.protein\nshow surface, polymer.protein\nset transparency, 0.35\nshow sticks, resn '+ligand+'\ncolor '+('teal' if ligand=='DPH' else 'orange')+', resn '+ligand+'\nselect contacts, byres (polymer.protein within 4.5 of resn '+ligand+')\nshow sticks, contacts\nzoom all\nbg_color white\nset ray_opaque_background, off\n')
 print(json.dumps(summary,indent=2),flush=True)
if __name__ == '__main__':analyze('PROTEIN_DPH_R3','DPH')
