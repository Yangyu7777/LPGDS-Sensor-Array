from pathlib import Path
import sys,json,csv,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageOps,ImageDraw
from openmm import app
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT;FIG=OUT/'figures';DATA=OUT/'analysis';FIG.mkdir(exist_ok=True);DATA.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'scripts'))
from md_analysis_utils import read_dcd,kabsch
SRC=ROOT/'analysis';RAW=ROOT/'raw'
def read(p):return np.genfromtxt(p,delimiter=',',names=True,encoding='utf-8')
m=read(SRC/'metrics.csv');s=read(SRC/'DPH_SASA.csv');w=read(SRC/'DPH_solvent.csv');r=read(SRC/'protein_RMSF.csv');c=np.genfromtxt(SRC/'residue_contacts.csv',delimiter=',',names=True,dtype=None,encoding='utf-8');t=m['production_time_ns']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
blue='#176781';gold='#c07835';gray='#6b7078'
figs=[]
def save(fig,name):
 for a in fig.axes:a.grid(False)
 for ext in ['png','pdf','svg']:fig.savefig(FIG/f'{name}.{ext}',dpi=400,bbox_inches='tight',facecolor='white')
 figs.append(name);plt.close(fig)
def axes(n,title):
 f,a=plt.subplots(n,1,figsize=(7.2,2.3*n),layout='constrained',squeeze=False);a=a[:,0]
 for x in a:x.set(xlabel='Production time (ns)',xlim=(0,5))
 f.suptitle(title);return f,a
f,a=axes(3,'FA R2 | Pairwise proximity')
for ax,key,label,color in zip(a,['DPH_protein_min_heavy_A','competitor_protein_min_heavy_A','DPH_competitor_min_heavy_A'],['DPH–WT','C18:2–WT','DPH–C18:2'],[blue,gold,gray]):
 ax.plot(t,m[key],color=color,lw=.9);ax.axhline(4.5,color='gray',ls='--',lw=.8,label='Contact cutoff: 4.5 Å');ax.set_ylabel(label+' distance (Å)');ax.legend(frameon=False,fontsize=8)
save(f,'Fig01_pairwise_distances')
# Supplement existing observables from real solute coordinates, without new MD.
p=app.PDBFile(str(RAW/'solute-topology.pdb'));at=list(p.topology.atoms());x,b=read_dcd(RAW/'solute-trajectory.dcd');x=x[500:];b=b[500:];assert len(x)==5000
ph=np.array([a.index for a in at if a.residue.name not in ['DPH','FAM'] and a.element.symbol!='H']);ca=np.array([a.index for a in at if a.name=='CA' and a.residue.name not in ['DPH','FAM']]);fit=ca[:170];dh=np.array([a.index for a in at if a.residue.name=='DPH' and a.element.symbol!='H']);fh=np.array([a.index for a in at if a.residue.name=='FAM' and a.element.symbol!='H']);mass=np.array([at[i].element.mass.value_in_unit(__import__('openmm').unit.dalton) for i in ph])
ref=np.array(app.PDBFile(str(ROOT/'snapshots/aligned_snapshot_0000ps.pdb')).positions.value_in_unit(__import__('openmm').unit.angstrom));vals=[]
for q,box in zip(x,b):
 rot,mc,tc=kabsch(q[fit],ref[fit]);aligned=(q-mc)@rot+tc;center=np.average(q[ph],axis=0,weights=mass);rg=np.sqrt(np.average(np.sum((q[ph]-center)**2,axis=1),weights=mass));v=[rg]
 for sel in [dh,fh]:
  mol=q[sel]-np.round((q[sel].mean(axis=0)-mc)/box)*box;fa=(mol-mc)@rot+tc;rr,mm,tt=kabsch(mol,ref[sel]);v += [np.sqrt(np.mean(np.sum((fa-ref[sel])**2,axis=1))),np.sqrt(np.mean(np.sum(((mol-mm)@rr+tt-ref[sel])**2,axis=1))),np.linalg.norm(fa.mean(axis=0)-ref[sel].mean(axis=0))]
 vals.append(v)
v=np.array(vals);np.savetxt(DATA/'additional_motion.csv',np.c_[t,v],delimiter=',',header='production_time_ns,protein_heavy_mass_weighted_Rg_A,DPH_protein_fit_RMSD_A,DPH_self_fit_RMSD_A,DPH_centroid_displacement_A,FA_protein_fit_RMSD_A,FA_self_fit_RMSD_A,FA_centroid_displacement_A',comments='')
f,a=axes(2,'FA R2 | Protein stability')
a[0].plot(t,m['protein_core_CA_RMSD_source_endpoint_A'],color=blue);a[0].set_ylabel('Core Cα RMSD (Å)');a[1].plot(t,v[:,0],color=blue);a[1].set_ylabel('Protein heavy-atom Rg (Å)');save(f,'Fig02_protein_RMSD_and_Rg')
f,a=plt.subplots(figsize=(8,3.5),layout='constrained');a.plot(np.arange(1,len(r)+1),r['CA_RMSF_A'],color=blue);a.axvspan(170.5,176.5,color='#eeeeee',label='His6 tag');a.set(xlabel='Input residue number',ylabel='Cα RMSF (Å)',title='FA R2 | Protein flexibility');a.legend(frameon=False);save(f,'Fig03_protein_RMSF')
f,a=axes(3,'FA R2 | Ligand motion relative to the source endpoint')
for ax,y1,y2,label in [(a[0],v[:,1],v[:,4],'Protein-aligned RMSD (Å)'),(a[1],v[:,2],v[:,5],'Self-aligned RMSD (Å)'),(a[2],v[:,3],v[:,6],'Centroid displacement (Å)')]:
 ax.plot(t,y1,color=blue,label='DPH',lw=.8);ax.plot(t,y2,color=gold,label='C18:2',lw=.8);ax.set_ylabel(label);ax.legend(frameon=False)
save(f,'Fig04_ligand_motion')
f,a=axes(3,'FA R2 | DPH surface exposure and solvation')
a[0].plot(s['production_time_ns'],s['DPH_isolated_SASA_A2'],ls='--',color=gray,label='Isolated DPH, same conformation');a[0].plot(s['production_time_ns'],s['DPH_protein_only_SASA_A2'],color=blue,label='WT occlusion only');a[0].plot(s['production_time_ns'],s['DPH_protein_plus_competitor_SASA_A2'],color=gold,label='WT + C18:2 occlusion');a[0].set_ylabel('DPH SASA (Å²)');a[0].legend(frameon=False,fontsize=8)
a[1].plot(s['production_time_ns'],100*s['DPH_protein_occluded_fraction'],color=blue,label='WT only');a[1].plot(s['production_time_ns'],100*s['DPH_total_occluded_fraction'],color=gold,label='WT + C18:2');a[1].set(ylabel='Occluded DPH SASA (%)',ylim=(0,100));a[1].legend(frameon=False)
a[2].plot(w['production_time_ns'],w['water_O_within5A'],color=blue,lw=.8);a[2].set_ylabel('Water O within 5 Å of DPH');save(f,'Fig05_DPH_exposure_and_hydration')
f,a=axes(2,'FA R2 | Local contacts and DMSO')
a[0].plot(t,m['DPH_contact_residues'],color=blue,label='DPH');a[0].plot(t,m['competitor_contact_residues'],color=gold,label='C18:2');a[0].set_ylabel('Contacting WT residues');a[0].legend(frameon=False);a[1].plot(w['production_time_ns'],w['DMSO_within4_5A'],color=blue);a[1].set_ylabel('DMSO within 4.5 Å of DPH');save(f,'Fig06_contacts_and_DMSO')
z=np.load(SRC/'contacts.npz');labels=z['residue_labels'];f,a=plt.subplots(1,2,figsize=(10,4.7),layout='constrained')
for ax,key,title,cmap in zip(a,['DPH_contact','competitor_contact'],['DPH–WT residue contacts','C18:2–WT residue contacts'],['Blues','YlOrBr']):
 mat=z[key];order=np.argsort(-mat.mean(axis=0),kind='stable')[:12];bins=mat.reshape(100,50,-1).mean(axis=1).T;im=ax.imshow(100*bins[order],aspect='auto',origin='upper',extent=[0,5,11.5,-.5],vmin=0,vmax=100,cmap=cmap);ax.set(yticks=range(12),yticklabels=labels[order],xlabel='Production time (ns)',title=title);f.colorbar(im,ax=ax,label='Contact occupancy (%)',shrink=.8)
save(f,'Fig07_residue_contact_heatmaps')
f,a=plt.subplots(1,2,figsize=(9,4.5),layout='constrained')
for ax,key,title,color in zip(a,['DPH_contact_fraction','competitor_contact_fraction'],['DPH–WT','C18:2–WT'],[blue,gold]):
 order=np.argsort(-c[key],kind='stable')[:12][::-1];ax.barh(c['input_residue'][order],100*c[key][order],color=color);ax.set(xlabel='Contact occupancy (%)',xlim=(0,100),title=title)
save(f,'Fig08_residue_contact_occupancy')
th=read(RAW/'thermodynamics.csv');sel=th['Time_ps']>502 if 'Time_ps' in th.dtype.names else None
# OpenMM quoted reporter header normalizes punctuation via genfromtxt.
print('THERMO_COLUMNS',th.dtype.names)
keys=th.dtype.names;timekey=next(k for k in keys if k.startswith('Time'));temp=next(k for k in keys if k.startswith('Temperature'));den=next(k for k in keys if k.startswith('Density'));mask=th[timekey]>502
f,a=axes(2,'FA R2 | Production temperature and density');a[0].plot((th[timekey][mask]-502)/1000,th[temp][mask],color=blue,lw=.7);a[0].axhline(298.15,color=gray,ls='--');a[0].set_ylabel('Temperature (K)');a[1].plot((th[timekey][mask]-502)/1000,th[den][mask],color=blue,lw=.7);a[1].set_ylabel('Density (g mL⁻¹)');save(f,'FigS01_temperature_and_density')
# Real, aligned coordinate snapshots; not a simulated or illustrated cavity surface.
f,a=plt.subplots(1,4,figsize=(12,3.6),layout='constrained',subplot_kw={'projection':'3d'})
for ax,ps,title in zip(a,[0,502,3002,5502],['Before preparation','Production 0 ns','Production 2.5 ns','Production 5 ns']):
 sp=app.PDBFile(str(ROOT/'snapshots'/f'aligned_snapshot_{ps:04d}ps.pdb'));q=np.array(sp.positions.value_in_unit(__import__('openmm').unit.angstrom));ax.plot(*q[ca[:170]].T,color='#aaaaaa',lw=1,alpha=.75)
 for sel,color,label in [(dh,blue,'DPH'),(fh,gold,'C18:2')]:
  for u,j in sp.topology.bonds():
   if u.index in sel and j.index in sel:ax.plot(*q[[u.index,j.index]].T,color=color,lw=2)
  ax.scatter(*q[sel].T,s=7,color=color)
 center=ref[ca[:170]].mean(axis=0)
 for func,mid in zip([ax.set_xlim,ax.set_ylim,ax.set_zlim],center):func(mid-26,mid+26)
 ax.view_init(elev=15,azim=-60);ax.set_box_aspect((1,1,1));ax.set_axis_off();ax.set_title(title,fontsize=10)
f.suptitle('FA R2 | Aligned coordinates: WT Cα trace (gray), DPH (blue), C18:2 (orange)');save(f,'FigS02_coordinate_snapshots')
summary=json.loads((SRC/'summary.json').read_text());summary.update(protein_core_CA_RMSD_mean_A=float(m['protein_core_CA_RMSD_source_endpoint_A'].mean()),protein_Rg_mean_A=float(v[:,0].mean()),protein_core_CA_RMSF_mean_A=float(r['CA_RMSF_A'][:170].mean()),DPH_mean_self_fit_RMSD_A=float(v[:,2].mean()),FA_mean_self_fit_RMSD_A=float(v[:,5].mean()),DPH_mean_DMSO_neighbors=float(w['DMSO_within4_5A'].mean()),DPH_mean_total_occluded_fraction=float(s['DPH_total_occluded_fraction'].mean()),mean_temperature_K=float(th[temp][mask].mean()),mean_density_g_ml=float(th[den][mask].mean()))
(DATA/'presentation_summary.json').write_text(json.dumps(summary,indent=2))
thumbs=[]
for i,n in enumerate(figs):
 im=Image.open(FIG/f'{n}.png').convert('RGB');im.thumbnail((520,440));tile=Image.new('RGB',(560,480),'white');tile.paste(im,((560-im.width)//2,24));ImageDraw.Draw(tile).text((15,460),n,fill='black');thumbs.append(tile)
canvas=Image.new('RGB',(1120,480*((len(thumbs)+1)//2)), '#eeeeee')
for i,im in enumerate(thumbs):canvas.paste(im,((i%2)*560,(i//2)*480))
canvas.save(FIG/'contact_sheet.jpg',quality=90)
print(json.dumps(summary,indent=2));print('FIGURES',figs)
