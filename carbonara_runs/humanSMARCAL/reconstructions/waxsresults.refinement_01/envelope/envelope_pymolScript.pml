load prot+solvlayer_0.pdb, input
run envelope_clipped.py
cmd.disable('envelope_origin')
cmd.hide("everything","input")

select water, resn sol+hoh
create waterObj, water


select rest, not (water or waterObj)

show sticks, waterObj
# set stick_radius, 0.5, waterObj (for newer PyMol)
set stick_radius, 0.15, waterObj
set stick_transparency, 0.0, waterObj

cmd.show("spheres"   ,"rest")
set sphere_scale, 0.85, symbol H
util.cba(22,"rest",_self=cmd)

cmd.zoom()
bg_color white
ray 480,360
png envfig_simple.png
