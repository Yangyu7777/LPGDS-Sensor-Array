from pathlib import Path
import csv,json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw
R=Path(__file__).resolve().parents[1];P=R.parent;F=R/'figures';A=R/'analysis'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
m=np.genfromtxt(A/'metrics.csv',delimiter=',',names=True);s=np.genfromtxt(A/'DPH_SASA.csv',delimiter=',',names=True);w=np.genfromtxt(A/'DPH_solvent.csv',delimiter=',',names=True);j=json.loads((A/'presentation_summary.json').read_text());sel=json.loads((R/'provenance/contact_snapshot_selection.json').read_text());t=m['production_time_ns'];joint=(m['DPH_protein_min_heavy_A']<4.5)&(m['competitor_protein_min_heavy_A']<4.5)&(m['DPH_competitor_min_heavy_A']<4.5)
def save(fig,n):
 for ext in ['png','pdf','svg']:fig.savefig(F/(n+'.'+ext),dpi=400,bbox_inches='tight',facecolor='white')
 plt.close(fig)
f,a=plt.subplots(2,1,figsize=(7.2,5),layout='constrained');tt=np.arange(100)*.05+.025
for key,color,label in [('DPH_protein_min_heavy_A','#176781','DPH–WT'),('competitor_protein_min_heavy_A','#c07835','C18:2–WT'),('DPH_competitor_min_heavy_A','#6b7078','DPH–C18:2')]:a[0].plot(tt,100*(m[key]<4.5).reshape(100,50).mean(axis=1),color=color,label=label,lw=1)
a[0].set(ylabel='Contact occupancy (%)',ylim=(-3,103));a[0].legend(frameon=False,fontsize=8)
a[1].fill_between(t,0,joint.astype(int),step='mid',color='#176781');a[1].set(ylabel='All three pairs contact',yticks=[0,1],yticklabels=['No','Yes'],ylim=(-.05,1.08));a[1].axvline(sel['selected_production_time_ns'],ls='--',color='#c07835',lw=1,label='Selected illustration');a[1].legend(frameon=False,fontsize=8)
for ax in a:ax.set(xlabel='Production time (ns)',xlim=(0,5));ax.grid(False)
f.suptitle('WT–DPH–(9Z,12Z)-C18:2 | Simultaneous proximity');save(f,'Fig09_simultaneous_three_pair_contacts')
# Remove white margins using one common crop for chronological views; no coordinate changes.
paths=[F/'structure_views'/f'{n}.png' for n in ['initial','production_0ns','production_2_5ns','production_5ns','selected_contact']]
bounds=[]
for path in paths:
 im=np.asarray(Image.open(path).convert('RGB'));ys,xs=np.where((im<240).any(axis=2));bounds.append((xs.min(),ys.min(),xs.max()+1,ys.max()+1))
box=(max(0,min(b[0] for b in bounds)-20),max(0,min(b[1] for b in bounds)-20),min(1200,max(b[2] for b in bounds)+20),min(1000,max(b[3] for b in bounds)+20))
def molecular_image(name):
 im=Image.open(F/'structure_views'/f'{name}.png').convert('RGB')
 if name!='selected_contact_closeup':return im.crop(box)
 ar=np.asarray(im);ys,xs=np.where((ar<240).any(axis=2));return im.crop((max(0,xs.min()-20),max(0,ys.min()-20),min(im.width,xs.max()+21),min(im.height,ys.max()+21)))
f,axes=plt.subplots(1,4,figsize=(12,3.3),layout='constrained')
for ax,name,title in zip(axes,['initial','production_0ns','production_2_5ns','production_5ns'],['Before preparation','Production 0 ns','Production 2.5 ns','Production 5 ns']):ax.imshow(molecular_image(name));ax.axis('off');ax.set_title(title,fontsize=9)
f.suptitle('WT–DPH–(9Z,12Z)-C18:2 | WT gray, DPH blue, linoleate orange');save(f,'FigS03_molecular_snapshots')
f,axes=plt.subplots(1,2,figsize=(8,3.7),layout='constrained')
for ax,name,title in zip(axes,['selected_contact','selected_contact_closeup'],['Whole protein','Contact-region close-up']):ax.imshow(molecular_image(name));ax.axis('off');ax.set_title(title)
f.suptitle(f"Selected simultaneous-contact illustration | Production {sel['selected_production_time_ns']:.3f} ns");save(f,'Fig10_selected_contact_structure')
if (P/'runs/FA_R1/status.json').exists():
 rows=[]
 for rep in [1,2,3]:
  st=json.loads((P/f'runs/FA_R{rep}/status.json').read_text());f=P/f'exports/FA_R{rep}/summary.json';row={'system':'cis,cis linoleate','original_run':f'FA_R{rep}','status':st['status'],'total_time_ps':st.get('time_ps'),'DPH_contact_fraction':None,'FA_contact_fraction':None,'DPH_FA_contact_fraction':None,'joint_contact_fraction':None,'WT_occluded_DPH_fraction':None,'water_O_mean':None}
  if f.exists():
   z=json.loads(f.read_text());row.update(DPH_contact_fraction=z['DPH_contact_fraction'],FA_contact_fraction=z['competitor_contact_fraction'],DPH_FA_contact_fraction=z['DPH_competitor_contact_fraction'],joint_contact_fraction=z['simultaneous_three_pair_contact_fraction'],WT_occluded_DPH_fraction=z['DPH_mean_protein_occluded_fraction'],water_O_mean=z['DPH_mean_water_neighbors'])
  rows.append(row)
 for rep in [1,2,3]:
  z=json.loads((P.parent/f'DPH_R3_Competition_20261020/exports/CONTROL_R{rep}/summary.json').read_text());rows.append({'system':'WT–DPH endpoint-continuation control','original_run':f'CONTROL_R{rep}','status':'complete','total_time_ps':5502,'DPH_contact_fraction':z['DPH_contact_fraction'],'FA_contact_fraction':None,'DPH_FA_contact_fraction':None,'joint_contact_fraction':None,'WT_occluded_DPH_fraction':z['DPH_mean_protein_occluded_fraction'],'water_O_mean':z['DPH_mean_water_neighbors']})
else:
 with (A/'repeat_and_control_context.csv').open() as h:rows=list(csv.DictReader(h))
with (A/'repeat_and_control_context.csv').open('w') as h:wr=csv.DictWriter(h,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
blocks=[]
for b in range(5):
 mm=(t>b)&(t<=b+1);ss=(s['production_time_ns']>b)&(s['production_time_ns']<=b+1);ww=(w['production_time_ns']>b)&(w['production_time_ns']<=b+1);blocks.append([b,b+1,int(mm.sum()),joint[mm].mean(),m['DPH_protein_min_heavy_A'][mm].mean(),m['competitor_protein_min_heavy_A'][mm].mean(),m['DPH_competitor_min_heavy_A'][mm].mean(),s['DPH_protein_occluded_fraction'][ss].mean(),s['DPH_total_occluded_fraction'][ss].mean(),w['water_O_within5A'][ww].mean()])
with (A/'one_ns_block_descriptive_means.csv').open('w') as h:wr=csv.writer(h);wr.writerow(['start_ns','end_ns','solute_frames','joint_contact_fraction','DPH_WT_distance_mean_A','FA_WT_distance_mean_A','DPH_FA_distance_mean_A','WT_occluded_DPH_mean','WT_FA_occluded_DPH_mean','water_O_mean']);wr.writerows(blocks)
figs=sorted(F.glob('Fig*.png'));canvas=Image.new('RGB',(1120,480*((len(figs)+1)//2)),'#eeeeee')
for i,path in enumerate(figs):
 im=Image.open(path).convert('RGB');im.thumbnail((530,440));tile=Image.new('RGB',(560,480),'white');tile.paste(im,((560-im.width)//2,15));ImageDraw.Draw(tile).text((10,458),path.stem,fill='black');canvas.paste(tile,((i%2)*560,(i//2)*480))
canvas.save(F/'contact_sheet.jpg',quality=90)
print('Figure sets:',len(figs));print('Repeat status:',[(r['original_run'],r['status']) for r in rows])
