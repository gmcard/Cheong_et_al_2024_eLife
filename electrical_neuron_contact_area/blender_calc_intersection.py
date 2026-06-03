import bpy
import os
import bpy
import bmesh

write_path=open('C:/path/to/manc_meshes/mn_elec_in_intersection_area.csv','w+')
write_path.write('mn,in,area,\n')

mn_path = 'C:/path/to/manc_meshes/mn_fattened/' 
elec_in_path = 'C:/path/to/manc_meshes/elec_in_fattened/'
intersection_path = 'C:/path/to/manc_meshes/mn_in_intersection/'

for obj in bpy.data.objects:
    bpy.data.objects.remove(obj, do_unlink=True)
            
for root, dirs, files in os.walk(mn_path):
    for f in files:
        if f.endswith('.obj') :
            for root, dirs, files_2 in os.walk(elec_in_path):
                for g in files_2:
                    if g.endswith('.obj') :
                        bpy.ops.import_scene.obj(filepath=os.path.join(mn_path, f))
                        bpy.ops.import_scene.obj(filepath=os.path.join(elec_in_path, g))
                        bpy.context.view_layer.objects.active = bpy.data.objects[0]
                        #bpy.ops.object.mode_set(mode='EDIT')
                        #bpy.ops.mesh.intersect_boolean(operation='INTERSECT', threshold = 1, solver='FAST')

                        bpy.ops.object.modifier_add(type='BOOLEAN')
                        bpy.context.object.modifiers["Boolean"].operation = 'INTERSECT'
                        bpy.context.object.modifiers["Boolean"].solver = 'EXACT'
                        #bpy.context.object.modifiers["Boolean"].double_threshold = 0.0001
                        bpy.context.object.modifiers["Boolean"].use_hole_tolerant = True
                        bpy.context.object.modifiers["Boolean"].object = bpy.data.objects[1]
                        bpy.ops.object.modifier_apply(modifier="Boolean")

                        #select the modified object and calculate area
                        bpy.context.view_layer.objects.active = bpy.data.objects[0]
                        obj = bpy.context.active_object
                        bm = bmesh.new()
                        bm.from_mesh(obj.data)
                        area = sum(h.calc_area() for h in bm.faces)
                        bm.free()
                        
                        #if area is positive, save intersection mesh so we can look at it
                        if area > 0:
                            #bpy.context.view_layer.objects.active = bpy.data.objects[1]
                            #bpy.data.objects[1].select = True
                            #bpy.data.objects[0].select = False
                            #bpy.ops.object.delete()
                            bpy.data.objects.remove(bpy.data.objects[1], do_unlink=True)
                            bpy.ops.export_scene.obj(filepath=os.path.join(intersection_path, f.split('.')[0] + '-' + g.split('.')[0] + '.obj'), use_materials=False)
                        
                        #write area info
                        write_path.write(f.split('.')[0] + ',' + g.split('.')[0] + ',' + str(area) + ',\n')
                        #bpy.ops.object.select_all(action='SELECT')
                        #bpy.ops.object.delete()
                        for obj in bpy.data.objects:
                            bpy.data.objects.remove(obj, do_unlink=True)
                        for block in bpy.data.meshes:
                            if block.users == 0:
                                bpy.data.meshes.remove(block)
        #bpy.ops.wm.read_factory_settings(use_empty=True)
write_path.close()