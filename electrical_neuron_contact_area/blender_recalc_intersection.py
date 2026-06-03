import bpy
import os
import bpy
import bmesh

intersection_path = 'C:path/to/manc_meshes/mn_in_intersection/'

write_path=open('C:/path/to/manc_meshes/mn_elec_in_intersection_area_recalc.csv','w')
write_path.write('mn,in,area,\n')

for root, dirs, files in os.walk(intersection_path):
    for f in files:
        if f.endswith('.obj') :
			bpy.ops.import_scene.obj(filepath=os.path.join(intersection_path, f))
			
			#select the modified object and calculate area
			obj = bpy.data.objects[0]
			bm = bmesh.new()
			bm.from_mesh(obj.data)
			area = sum(h.calc_area() for h in bm.faces)
			bm.free()
			
			f = f.split('.')[0]
			
			#write area info
			write_path.write(f.split('-')[0] + ',' + f.split('-')[1] + ',' + str(area) + ',\n')

			for obj in bpy.data.objects:
				bpy.data.objects.remove(obj, do_unlink=True)
                
write_path.close()

