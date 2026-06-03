import bpy
import os

path = 'C:/path/to/manc_meshes/mn/'
write_path = 'C:/path/to/manc_meshes/mn_fattened/'

for root, dirs, files in os.walk(path):
    for f in files:
        if f.endswith('.obj') :
            original_file = os.path.join(path, f)
            new_file = os.path.join(write_path, f)
            
            #bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.select_all(action='SELECT')
            bpy.ops.object.delete()

            bpy.ops.import_scene.obj(filepath=original_file)
            #need to set active object
            bpy.context.view_layer.objects.active = bpy.context.window.scene.objects[0]
            bpy.ops.object.mode_set(mode='EDIT', toggle=False)
            bpy.ops.mesh.select_all(action = 'SELECT')
            bpy.ops.transform.shrink_fatten(value = 0.1) #fatten by 100 nm
            
            bpy.ops.object.mode_set(mode='OBJECT', toggle=False)
            bpy.ops.object.select_all(action='SELECT')
            bpy.ops.export_scene.obj(filepath=new_file, use_materials=False)
            
            
            