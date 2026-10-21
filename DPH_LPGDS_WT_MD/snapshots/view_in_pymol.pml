reinitialize
load aligned_5502ps.pdb, final
hide everything
show cartoon, polymer.protein
color gray70, polymer.protein
show surface, polymer.protein
set transparency, 0.35
show sticks, resn DPH
color teal, resn DPH
select contacts, byres (polymer.protein within 4.5 of resn DPH)
show sticks, contacts
zoom all
bg_color white
set ray_opaque_background, off
