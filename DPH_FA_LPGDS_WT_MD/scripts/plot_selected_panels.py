from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
R=Path(__file__).resolve().parents[1];O=R/'figures/Selected_three_panels';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':False,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
def save(f,n):
 for e in ['png','pdf','svg']:f.savefig(O/f'{n}.{e}',dpi=600,bbox_inches='tight',facecolor='white')
 plt.close(f)
# Four chronological states, identical camera and identical white-margin crop.
names=['initial','production_0ns','production_2_5ns','production_5ns'];imgs=[Image.open(R/'figures/structure_views'/f'{n}.png').convert('RGB') for n in names];
for k,im in enumerate(imgs):
 ar=np.asarray(im).astype(float)/255
 import matplotlib.colors as mc
 hsv=mc.rgb_to_hsv(ar);mask=(hsv[:,:,0]>.45)&(hsv[:,:,0]<.66)&(hsv[:,:,1]>.22)
 fa=(hsv[:,:,0]>.035)&(hsv[:,:,0]<.14)&(hsv[:,:,1]>.22)
 hsv[:,:,0][mask]=.085
 hsv[:,:,0][fa]=.55
 imgs[k]=Image.fromarray(np.uint8(np.clip(mc.hsv_to_rgb(hsv),0,1)*255))
bs=[]
for im in imgs:
 a=np.array(im);y,x=np.where((a<240).any(axis=2));bs.append((x.min(),y.min(),x.max()+1,y.max()+1))
b=(max(0,min(x[0] for x in bs)-20),max(0,min(x[1] for x in bs)-20),min(1200,max(x[2] for x in bs)+20),min(1000,max(x[3] for x in bs)+20))
f,ax=plt.subplots(1,4,figsize=(10,3),layout='constrained')
for a,im,title in zip(ax,imgs,['Before preparation','0 ns','2.5 ns','5 ns']):a.imshow(im.crop(b));a.axis('off');a.set_title(title,fontsize=10)
f.text(.5,.01,'WT: gray   |   DPH: orange   |   cis-Linoleate: blue (O: red)',ha='center',fontsize=9);save(f,'01_Structural_evolution')
s=np.genfromtxt(R/'analysis/DPH_SASA.csv',delimiter=',',names=True);w=np.genfromtxt(R/'analysis/DPH_solvent.csv',delimiter=',',names=True)
f,ax=plt.subplots(3,1,figsize=(4.4,6.2),layout='constrained')
ax[0].plot(s['production_time_ns'],s['DPH_protein_only_SASA_A2'],color='#777777',lw=1.1,ls='--',label='DPH + protein');ax[0].plot(s['production_time_ns'],s['DPH_protein_plus_competitor_SASA_A2'],color='#c07835',lw=1.2,label='DPH + protein + linoleate');ax[0].set_ylabel('DPH SASA (Å²)')
ax[1].plot(s['production_time_ns'],s['DPH_protein_occluded_fraction']*100,color='#777777',lw=1.1,ls='--');ax[1].plot(s['production_time_ns'],s['DPH_total_occluded_fraction']*100,color='#c07835',lw=1.2);ax[1].set(ylabel='Occluded DPH SASA (%)',ylim=(0,100))
ax[2].plot(w['production_time_ns'],w['water_O_within5A'],color='#c07835',lw=.9);ax[2].set_ylabel('Water O within 5 Å of DPH')
for a in ax:a.set(xlim=(0,5),xlabel='Production time (ns)');a.grid(False)
save(f,'02_DPH_exposure_and_hydration')
z=np.load(R/'analysis/contacts.npz');mat=z['DPH_contact'];order=np.argsort(-mat.mean(axis=0),kind='stable')[:12];bins=mat.reshape(100,50,-1).mean(axis=1).T
f,a=plt.subplots(figsize=(4.5,6.2),layout='constrained');im=a.imshow(100*bins[order],aspect='auto',extent=[0,5,11.5,-.5],vmin=0,vmax=100,cmap='Oranges');a.set(yticks=range(12),yticklabels=z['residue_labels'][order],xlabel='Production time (ns)',title='DPH–WT contacts (with cis-linoleate)');f.colorbar(im,ax=a,location='top',label='Contact occupancy (%)',shrink=.9,pad=.03);save(f,'03_DPH_WT_contact_heatmap')
(O/'README.md').write_text('''# Selected R2 panels

Source: new cis,cis-linoleate batch, CAS60-33-3, original run FA_R2. PNG600dpi/PDF/SVG; no grid or other-run comparison.

1. Structural evolution: actual four chronological snapshots with a fixed camera/common crop. The pre-preparation frame precedes 502ps preparation; 0/2.5/5ns are production times. His6 hidden only in rendering. Molecular renders are raster images embedded in PDF/SVG.
2. DPH exposure/hydration: SASA in the complete WT+cis-linoleate environment, and combined WT+linoleate occluded fraction relative to isolated DPH of the same conformation. No isolated or protein-only comparison curve is drawn. Water counts use O atoms within5Å of DPH heavy atoms. SASA50ps; water10ps.
3. Heatmap: DPH–WT residue contacts within the ternary trajectory (not a separate binary control and not DPH–linoleate). Any intermolecular heavy-atom periodic distance<4.5Å; each cell is occupancy in50ps.

The reference image was from the older binary WT–DPH dataset. These panels use the selected current ternary R2 and therefore need not have the same curve shapes. No simulation data were changed.
''')
print(O)
