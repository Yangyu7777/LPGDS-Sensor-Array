from pathlib import Path
import json
import pymol2
R=Path(__file__).resolve().parents[1];out=R/'figures/structure_views';out.mkdir(exist_ok=True);sel=json.loads((R/'provenance/contact_snapshot_selection.json').read_text());files=[('initial','aligned_snapshot_0000ps.pdb'),('production_0ns','aligned_snapshot_0502ps.pdb'),('production_2_5ns','aligned_snapshot_3002ps.pdb'),('production_5ns','aligned_snapshot_5502ps.pdb'),('selected_contact',sel['snapshot'])]
with pymol2.PyMOL() as pm:
 c=pm.cmd;c.set('max_threads',2);c.bg_color('white');c.set('orthoscopic',1);c.set('ray_shadows',0);c.set('antialias',2);c.set('specular',.15);c.set('ambient',.5);c.set('cartoon_fancy_helices',1);c.set('stick_radius',.18);view=None
 for name,f in files:
  c.delete('all');c.load(str(R/'snapshots'/f),'system');c.hide('everything');c.show('cartoon','polymer.protein and not resi 171-176');c.color('gray75','polymer.protein');c.show('sticks','resn DPH+FAM and not elem H');c.set_color('dph_color',[.09,.40,.51]);c.set_color('fa_color',[.75,.47,.21]);c.color('dph_color','resn DPH');c.color('fa_color','resn FAM');c.color('red','resn FAM and elem O')
  if view is None:c.orient('polymer.protein and not resi 171-176');c.zoom('all',4);c.turn('y',30);view=c.get_view()
  c.set_view(view);c.png(str(out/(name+'.png')),width=1200,height=1000,dpi=400,ray=1);print('Rendered',name,flush=True)
 c.set('cartoon_transparency',.55);c.zoom('resn DPH+FAM',6);c.png(str(out/'selected_contact_closeup.png'),width=1200,height=1000,dpi=400,ray=1)
 c.save(str(out/'WT_DPH_linoleate_ZZ.pse'))
(R/'provenance/pymol_render.json').write_text(json.dumps({'camera':list(view),'same_camera_all_frames':True,'protein_render':'cartoon; His6 residues171-176 hidden only for visualization','ligands':'actual heavy-atom sticks; DPH blue, cis-linoleate orange with red oxygen','artificial_ligand_or_protein_movement':False},indent=2)+'\n')
