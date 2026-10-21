"""Export conditional competition observables; no outcome is selected or forced."""
from pathlib import Path
import os,sys,json,csv,hashlib
import numpy as np
from scipy.spatial import cKDTree
from openmm import app,unit
from md_analysis_utils import read_dcd,kabsch
root=Path(__file__).resolve().parents[1];mode='FA';rep=2;name='WT–DPH–C18:2';raw=root/'raw';out=root/'analysis';out.mkdir(exist_ok=True);source=root/'inputs'
assert json.loads((raw/'status.json').read_text())['status']=='complete'
p=app.PDBFile(str(raw/'solute-topology.pdb'));at=list(p.topology.atoms());xyz,boxes=read_dcd(raw/'solute-trajectory.dcd');time=np.arange(3,5503);assert len(xyz)==len(time);prod=time>502
ip=app.PDBFile(str(source/'solvated.pdb'));ia=list(ip.topology.atoms());keep=[a.index for a in ia if a.residue.name not in {'HOH','DMS','PH1','PH2','NA','K','CL'}];ref=np.load(source/'coordinates_nm.npy')[keep]*10
protein=np.array([a.index for a in at if a.residue.name not in ['DPH','FAM','PFO']]);ph=np.array([i for i in protein if at[i].element.symbol!='H']);dpall=np.array([a.index for a in at if a.residue.name=='DPH']);dh=np.array([i for i in dpall if at[i].element.symbol!='H']);comp=np.array([a.index for a in at if a.residue.name in ['FAM','PFO']],int);ch=np.array([i for i in comp if at[i].element.symbol!='H'],int);ca=np.array([i for i in protein if at[i].name=='CA']);fit=ca[:170];res=[r for r in p.topology.residues() if r.name not in ['DPH','FAM','PFO']];rmap={r.index:i for i,r in enumerate(res)};labels=[r.name+r.id for r in res];contact=np.zeros((len(xyz),len(res)),bool);cc=contact.copy();dc0=ref[dh].mean(axis=0);direction=np.array(json.loads((root/'provenance/preparation.json').read_text())['initial_outward_direction']);rows=[];aligned_ca=[];snap={0:ref.copy()};sasa_rows=[]
rad=np.array([{'H':1.2,'C':1.7,'N':1.55,'O':1.52,'S':1.8,'F':1.47}[a.element.symbol]+1.4 for a in at]);nn=960;idx=np.arange(nn);zz=1-2*(idx+.5)/nn;phi=idx*np.pi*(3-np.sqrt(5));rr=np.sqrt(1-zz*zz);directions=np.c_[rr*np.cos(phi),rr*np.sin(phi),zz]
def dph_sasa(x,b):
 iso=byprotein=total=0.
 for j,atom in enumerate(dpall):
  points=x[atom]+rad[atom]*directions;delta=points[:,None,:]-x[dpall][None,:,:];delta-=np.round(delta/b)*b;blocked=np.sum(delta**2,axis=2)<rad[dpall][None,:]**2-1e-8;blocked[:,j]=False;free=~blocked.any(axis=1);occ=[]
  for ob in [protein,comp]:
   if not len(ob):occ.append(np.zeros(nn,bool));continue
   delta=x[ob]-x[atom];delta-=np.round(delta/b)*b;near=np.flatnonzero(np.linalg.norm(delta,axis=1)<rad[atom]+rad[ob]);q=points[:,None,:]-(x[atom]+delta[near])[None,:,:];occ.append((np.sum(q*q,axis=2)<rad[ob[near]][None,:]**2-1e-8).any(axis=1))
  weight=4*np.pi*rad[atom]**2/nn;iso+=free.sum()*weight;byprotein+=(free&~occ[0]).sum()*weight;total+=(free&~occ[0]&~occ[1]).sum()*weight
 return [iso,byprotein,total,1-byprotein/iso,1-total/iso]
for k,(x,b,t) in enumerate(zip(xyz,boxes,time)):
 x=x.astype(float);tree=cKDTree(x[ph]%b,boxsize=b);dpmin=float(tree.query(x[dh]%b)[0].min());rot,mc,tc=kabsch(x[fit],ref[fit]);f=(x-mc)@rot+tc
 for sel in [dpall,comp]:
  if len(sel):shift=np.round((x[sel].mean(axis=0)-mc)/b)*b;f[sel]=(x[sel]-shift-mc)@rot+tc
 def touched(selection):
  if not len(selection):return set()
  near=cKDTree(x[selection]%b,boxsize=b).query(x[ph]%b)[0]<4.5;return {rmap[at[ph[j]].residue.index] for j in np.flatnonzero(near)}
 hits=touched(dh);contact[k,list(hits)]=True;hc=touched(ch);cc[k,list(hc)]=True
 if len(ch):cm=float(tree.query(x[ch]%b)[0].min());pair=float(cKDTree(x[dh]%b,boxsize=b).query(x[ch]%b)[0].min());occmin=float(cKDTree(ref[dh]).query(f[ch])[0].min());ccent=float(np.linalg.norm(f[ch].mean(axis=0)-dc0))
 else:cm=pair=occmin=ccent=float('nan')
 dcen=f[dh].mean(axis=0);prms=float(np.sqrt(np.mean(np.sum((f[fit]-ref[fit])**2,axis=1))));lrms=float(np.sqrt(np.mean(np.sum((f[dh]-ref[dh])**2,axis=1))))
 rows.append([t,(t-502)/1000,dpmin,len(hits),lrms,np.linalg.norm(dcen-dc0),np.dot(dcen-dc0,direction),prms,cm,pair,occmin,ccent,len(hc)]);aligned_ca.append(f[ca])
 if t in [502,3002,5502]:snap[int(t)]=f.copy()
 if t>502 and t%50==0:sasa_rows.append([t,(t-502)/1000,*dph_sasa(x,b)])
v=np.asarray(rows);ss=np.asarray(sasa_rows);ac=np.asarray(aligned_ca)[prod];rf=np.sqrt(np.mean(np.sum((ac-ac.mean(axis=0))**2,axis=2),axis=0))
def write(n,h,rows):
 with (out/n).open('w') as file:w=csv.writer(file);w.writerow(h);w.writerows(rows)
write('metrics.csv',['total_time_ps','production_time_ns','DPH_protein_min_heavy_A','DPH_contact_residues','DPH_protein_aligned_RMSD_source_endpoint_A','DPH_centroid_displacement_A','DPH_outward_projection_A','protein_core_CA_RMSD_source_endpoint_A','competitor_protein_min_heavy_A','DPH_competitor_min_heavy_A','competitor_to_original_DPH_atoms_min_A','competitor_centroid_to_original_DPH_centroid_A','competitor_contact_residues'],v[prod]);write('DPH_SASA.csv',['total_time_ps','production_time_ns','DPH_isolated_SASA_A2','DPH_protein_only_SASA_A2','DPH_protein_plus_competitor_SASA_A2','DPH_protein_occluded_fraction','DPH_total_occluded_fraction'],ss)
write('protein_RMSF.csv',['input_residue','CA_RMSF_A'],zip(labels,rf));write('residue_contacts.csv',['input_residue','DPH_contact_fraction','competitor_contact_fraction'],zip(labels,contact[prod].mean(axis=0),cc[prod].mean(axis=0)));np.savez_compressed(out/'contacts.npz',time_ps=time[prod],residue_labels=labels,DPH_contact=contact[prod],competitor_contact=cc[prod])
fp=app.PDBFile(str(raw/'full-topology.pdb'));fa=list(fp.topology.atoms());full,fb=read_dcd(raw/'full-trajectory.dcd');fd=[a.index for a in fa if a.residue.name=='DPH' and a.element.symbol!='H'];water=[a.index for a in fa if a.residue.name=='HOH' and a.element.symbol=='O'];dm=[a.index for a in fa if a.residue.name=='DMS' and a.element.symbol!='H'];sol=[]
for i,(x,b) in enumerate(zip(full,fb)):
 t=(i+1)*10
 if t<=502:continue
 tree=cKDTree(x[fd]%b,boxsize=b);wd=tree.query(x[water]%b)[0];dd=tree.query(x[dm]%b)[0];sol.append([t,(t-502)/1000,int((wd<5).sum()),len({fa[dm[j]].residue.index for j in np.flatnonzero(dd<4.5)})])
sol=np.asarray(sol);assert len(sol)==500;del full,fb
write('DPH_solvent.csv',['total_time_ps','production_time_ns','water_O_within5A','DMSO_within4_5A'],sol)
p.topology.setPeriodicBoxVectors(None)
for t,f in snap.items():
 with (root/'snapshots'/f'aligned_snapshot_{t:04d}ps.pdb').open('w') as file:app.PDBFile.writeFile(p.topology,f*unit.angstrom,file)
# Descriptive numerical results; no binding, displacement or cavity classification from one metric.
last=time>4502;j={'run':name,'branch':mode,'replica':rep,'production_frames':5000,'production_ns':5,'DPH_contact_fraction':float((v[prod,2]<4.5).mean()),'DPH_last_1ns_contact_fraction':float((v[last,2]<4.5).mean()),'DPH_final_min_protein_A':float(v[-1,2]),'DPH_final_centroid_displacement_A':float(v[-1,5]),'DPH_final_outward_projection_A':float(v[-1,6]),'DPH_mean_protein_occluded_fraction':float(ss[:,5].mean()),'DPH_last_saved_SASA_protein_occluded_fraction':float(ss[-1,5]),'SASA_last_time_ps':float(ss[-1,0]),'DPH_mean_water_neighbors':float(sol[:,2].mean()),'competitor_contact_fraction':float((v[prod,8]<4.5).mean()) if len(comp) else None,'competitor_final_min_protein_A':float(v[-1,8]) if len(comp) else None,'DPH_competitor_contact_fraction':float((v[prod,9]<4.5).mean()) if len(comp) else None,'scope':'Endpoint-conditional exploratory branch; inspect all repeats and control; no guaranteed replacement or closed-cavity classification','sha256':{n:hashlib.sha256((raw/n).read_bytes()).hexdigest() for n in ['solute-trajectory.dcd','full-trajectory.dcd']}}
(out/'summary.json').write_text(json.dumps(j,indent=2)+'\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':False,'pdf.fonttype':42,'svg.fonttype':'none'})
t=v[prod,1];col={'CONTROL':'#666666','FA':'#a75a32','PFOA':'#16657a'}[mode]
def save(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(out/(name+'.'+ext),dpi=350,bbox_inches='tight')
 plt.close(fig)
fig,axes=plt.subplots(3,1,figsize=(7.2,7.8),layout='constrained');axes[0].plot(t,v[prod,2],color=col);axes[0].axhline(4.5,color='gray',ls='--');axes[0].set_ylabel('DPH–protein min distance (Å)');axes[1].plot(t,v[prod,5],color=col,label='Centroid displacement');axes[1].plot(t,v[prod,6],color='gray',ls='--',label='Outward projection');axes[1].set_ylabel('DPH movement (Å)');axes[1].legend(frameon=False);axes[2].plot(t,v[prod,4],color=col,label='DPH, protein aligned');axes[2].plot(t,v[prod,7],color='gray',label='Protein core Cα');axes[2].set_ylabel('RMSD from source endpoint (Å)');axes[2].legend(frameon=False)
for a in axes:a.set(xlabel='Production time (ns)',xlim=(0,5));a.grid(False)
fig.suptitle(name+' | Resident DPH behavior');save(fig,'01_DPH_position_and_contacts')
fig,axes=plt.subplots(3,1,figsize=(7.2,7.8),layout='constrained');axes[0].plot(ss[:,1],ss[:,3],color=col,label='Protein occlusion only');axes[0].plot(ss[:,1],ss[:,4],color='gray',ls='--',label='Protein + competitor');axes[0].set_ylabel('DPH SASA (Å²)');axes[0].legend(frameon=False);axes[1].plot(ss[:,1],ss[:,5]*100,color=col);axes[1].set(ylabel='Protein-occluded DPH SASA (%)',ylim=(0,100));axes[2].plot(sol[:,1],sol[:,2],color=col);axes[2].set_ylabel('Water O within 5 Å of DPH')
for a in axes:a.set(xlabel='Production time (ns)',xlim=(0,5));a.grid(False)
fig.suptitle(name+' | Resident DPH exposure');save(fig,'02_DPH_exposure_and_water')
if len(comp):
 fig,axes=plt.subplots(3,1,figsize=(7.2,7.8),layout='constrained')
 for ax,y,label in [(axes[0],8,'Competitor–protein min distance (Å)'),(axes[1],9,'Competitor–DPH min distance (Å)'),(axes[2],10,'Competitor–initial DPH region (Å)')]:ax.plot(t,v[prod,y],color=col);ax.set(xlabel='Production time (ns)',ylabel=label,xlim=(0,5));ax.grid(False)
 fig.suptitle(name+' | Competitor approach (distance proxies)');save(fig,'03_competitor_approach')
fig,ax=plt.subplots(figsize=(8,3.8),layout='constrained');ax.plot(range(1,len(rf)+1),rf,color=col);ax.set(xlabel='Input residue number (His6 retained)',ylabel='Cα RMSF (Å)',title=name);save(fig,'04_protein_RMSF')
print(json.dumps(j,indent=2))
