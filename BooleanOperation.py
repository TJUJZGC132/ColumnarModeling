#-----------------------------------------------#
#General columnar grain boundary modeling method#
#-----------------------------------------------#
import rhinoscriptsyntax as rs
import scriptcontext as sc
import math
#----------------------------------------------#
#Define the boundary parameters of the specimen#
#----------------------------------------------#
b = 50
a = 50
h = 100
corner1 = [0, 0, 0]
corner2 = [a, 0, 0]
corner3 = [a, b, 0]
corner4 = [0, b, 0]
corner5 = [0, 0, h]
corner6 = [a, 0, h]
corner7 = [a, b, h]
corner8 = [0, b, h]
#--------------------------#
#Create a physical specimen#
#--------------------------#
box_corners = [corner1, corner2, corner3, corner4, corner5, corner6, corner7, corner8]
box = rs.AddBox(box_corners)
translation_vector = [50, 50, 50]
moved_box = rs.MoveObject(box, translation_vector)
#-------------------------------------#
#Create a Brazilian Splitting Specimen#
#-------------------------------------#
def create_cylinder(base_point, radius, height):
    cylinder = rs.AddCylinder(base_point, height, radius)
    return cylinder
#---------------#
#Create a sphere#
#---------------#
def create_closed_sphere(center, radius):
    sphere = rs.AddSphere(center, radius)
    return sphere
#-------------------------------#
#Input columnar grain boundaries#
#-------------------------------#
def import_nodes_from_txt(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        points = []

        for line in lines:
            coords = line.strip().split(' ')
            if len(coords) == 3:
                x, y, z = map(float, coords)
                point = rs.AddPoint(x, y, z)
                points.append(point)
            elif len(coords) == 2:
                x, y = map(float, coords)
                point = rs.AddPoint(x, y)
                points.append(point)
            else:
                print("Invalid line format: {}".format(line))

        return points

    except IOError:
        print("Error reading the file.")
#-------------------------#
#Generate pre-cut specimen#
#-------------------------#
def create_faces_from_points(points):
    faces = []  
    for i in range(0, len(points) - 3, 4):
        face = rs.AddSrfPt([points[i], points[i + 1], points[i + 2], points[i + 3]])
        if rs.IsBrep(face):
            if not is_duplicate_face(face, faces):
                faces.append(face)
            else:
                rs.DeleteObject(face)
        else:
            print("Invalid face geometry at index:", i)
    return faces
#-----------------------------------#
#Delete overlapping grain boundaries#
#-----------------------------------#
def is_duplicate_face(new_face, existing_faces):
    new_center = rs.SurfaceAreaCentroid(new_face)[0]
    
    for face in existing_faces:
        existing_center = rs.SurfaceAreaCentroid(face)[0]
        distance = rs.Distance(new_center, existing_center)
        
        if distance < 1e-6:
            return True

    return False
#---------------------#
#Define cutting method#
#---------------------#
def rhino_brep_split(brep_id, cutter_ids):
    rs.EnableRedraw(False)
    rs.UnselectAllObjects()
    
    rs.SelectObject(brep_id)
    cmd = "_-Split"
    for id in cutter_ids:
        cmd += " _SelID "
        cmd += id.ToString()
    cmd += " _Enter"
    rs.Command(cmd, False)
    results = rs.SelectedObjects()
    
    rs.UnselectAllObjects()
    rs.EnableRedraw(True)
    return results
#--------------------------#
#Define the cutting process#
#--------------------------#
def cut_faces_with_sphere(closed_sphere , face_list):
    cut_results = []

    for face in face_list:
        cut_result = rhino_brep_split(face, [closed_sphere ])
        cut_results.extend(cut_result)
    
    return cut_results
#----------------------------------#
#Delete irrelevant grain boundaries#
#----------------------------------#
def remove_intersecting_results_up(plant1, cut_results):
    remaining_results = [] 
    for cut_result in cut_results:
        intersection_curve = rs.IntersectBreps(cut_result, plant1)
        if intersection_curve:
            rs.DeleteObject(cut_result)
        else:
            remaining_result=cut_result
            remaining_results.append(remaining_result)
    return remaining_results
    
def remove_intersecting_results_down(plant2, remaining_results):
    final_results = [] 
    for remaining_result in remaining_results:
        intersection_curve = rs.IntersectBreps(remaining_result, plant2)
        if intersection_curve:
            rs.DeleteObject(remaining_result)
        else:
            final_result = remaining_result
            final_results.append(final_result)
    return final_results

def delete_all_lines():
    lines = rs.ObjectsByType(rs.filter.curve)
    for line in lines:
        if rs.IsCurve(line) and rs.IsLine(line):
            rs.DeleteObject(line)
#--------------------#
#Output triangle mesh#
#--------------------#
def export_face_vertices_to_file(final_results, output_file):
    try:
        with open(output_file, 'w') as file:
            for final_result in final_results:
                plane = final_result
                boundary_curves = rs.DuplicateEdgeCurves(plane)
                if boundary_curves:
                    closed_curve = rs.JoinCurves(boundary_curves, True)
                    if closed_curve:
                        vertices = rs.CurvePoints(closed_curve)
                        vertex_line = " ".join(["{} {} {}".format(vertex[0], vertex[1], vertex[2]) for vertex in vertices])
                        file.write(vertex_line + "\n") 
    except IOError:
        print("Error writing to the file.")

if __name__ == "__main__":
#--------------------------------------#
#Generate specimens of different shapes#
#--------------------------------------#
    sphere_center = (50, 50, 50)
    sphere_radius = 25
    closed_sphere = moved_box#Rectangular specimen
    #closed_sphere = create_closed_sphere(sphere_center, sphere_radius)#Spherical specimen
    #sphere_center = (1.0, 1.0, 1.0)
    #sphere_radius = 0.8
    #base_point = (50, 50, 25)  #Brazilian plate style parameters
    #cylinder_radius = 25  # Brazilian plate style parameters
    #cylinder_height = 25  # Brazilian plate style parameters
    #cutting_cylinder = create_cylinder(base_point, cylinder_radius, cylinder_height)#Brazilian split specimen
    #closed_sphere = cutting_cylinder#Brazilian split specimen
#-----------------#
#Define input file#
#-----------------#
    file_path = 'Your\\input\\file\\path.txt' 
    imported_points = import_nodes_from_txt(file_path)
    faces = create_faces_from_points(imported_points)
    rs.DeleteObjects(rs.ObjectsByType(rs.filter.point))
#-----------------------#
#Rotate the cut specimen#
#-----------------------#
    rotation_axis_x = [1, 0, 0]  # X
    rotation_axis_y = [0, 1, 0]  # Y
    rotation_axis_z = [0, 0, 1]  # Z
    angle_x = 0
    angle_y = 0
    angle_z = 0
    point_x = 75
    point_y = 75
    point_z = 75
    rs.RotateObject(faces, [point_x, point_y, point_z], angle_x, rotation_axis_x)
    rs.RotateObject(faces, [point_x, point_y, point_z], angle_y, rotation_axis_y)
    rs.RotateObject(faces, [point_x, point_y, point_z], angle_z, rotation_axis_z)
#--- ---------------------#
#Perform Boolean operations 
#-------------------------#
    cutting_sphere = closed_sphere  #”colsed_sphere” as the cutting surface
    face_list = faces #”faces” as the imported cut model
    cut_results = cut_faces_with_sphere(cutting_sphere , face_list)
#------------------------#
#Delete irrelevant planes#
#------------------------#
    h1=0
    h2=200
    h3=h2-h1
    corner1 = (0, 0, h1)
    corner2 = (h2, 0, h1)
    corner3 = (h2, h2, h1)
    corner4 = (0, h2, h1)
    corner5 = (0, 0, h3)
    corner6 = (h2, 0, h3)
    corner7 = (h2, h2, h3)
    corner8 = (0, h2, h3)
    plant1 = rs.AddSrfPt([corner1, corner2, corner3, corner4])
    plant2 = rs.AddSrfPt([corner5, corner6, corner7, corner8])
    rs.RotateObject(plant1, [point_x, point_y, point_z], angle_x, rotation_axis_x)
    rs.RotateObject(plant1, [point_x, point_y, point_z], angle_y, rotation_axis_y)
    rs.RotateObject(plant1, [point_x, point_y, point_z], angle_z, rotation_axis_z)
    
    rs.RotateObject(plant2, [point_x, point_y, point_z], angle_x, rotation_axis_x)
    rs.RotateObject(plant2, [point_x, point_y, point_z], angle_y, rotation_axis_y)
    rs.RotateObject(plant2, [point_x, point_y, point_z], angle_z, rotation_axis_z)
    
    
    plants = [plant1, plant2]
    remaining_results = remove_intersecting_results_up(plant1,cut_results)
    final_results=remove_intersecting_results_down(plant2,remaining_results)
    rs.DeleteObjects(plant1)
    rs.DeleteObjects(plant2)
    rs.DeleteObjects(cutting_sphere)
    delete_all_lines()
#-------------------------------#
#Return to the coordinate origin#
#-------------------------------#
    translation_vector_back = [-50, -50, -50]
    rs.MoveObject(final_results, translation_vector_back)
#------------------#
#Define output file#
#------------------#
    export_file_path = 'Your\\output\\file\\path.txt'
    export_face_vertices_to_file(final_results, export_file_path)
