import scipy.io as sio
from stl import mesh
from skimage import measure
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import subprocess
from scipy.spatial import ConvexHull
from solid import *
from solid.utils import *
import argparse
from shapely.geometry import Polygon, LineString
from scipy import ndimage

parser = argparse.ArgumentParser()
parser.add_argument("inputFile", type=str,
                    help="name of .mat input file containing 3D volume matrix")
#parser.add_argument("outputFile", type=str,
#                    help="name of .stl file which will contain final mould")
parser.add_argument("--show_convex_hull", help="displays the convex hull of tumor 2D projection")
args = parser.parse_args()

tumorVolumeNameInMat = 'threshold_smooth_box' # Design mechanism to select which matrix via prompt
tumorPointsNameInMat = 'labels_smooth_box' # Design mechanism to select which matrix via prompt
inputFilePath = args.inputFile
outputFilePath = "Patient"

# All units are given in millimeter
gapBetweenSlabs = 1.5 # Adapt to blade thickness
sectionPlaneDistances = 10 # Default in our case
sizeInXDirection = 210 # Printer limitation
sizeInYDirection = 210 # Printer limitation
sizeinZDirection = 250 # Printer limitation
basePlateThickness = 4 #
slabHeight = 60 # this
contactHoleSize = 10 # size of holes in print

sectionDirection = 1 # 0 = Along X?, 1 = Along Y?
tumorRotation = False # Not implemented, assume correct prior orientation

kidneyPointIndex = [2]
tumorPointIndex = [3]
vesselPointIndices = [6] # DODGY, REPLACE!!!
contactPointIndices = [7,8] # DODGY, REPLACE!!!
pointLabels = ['V','CK','CT']


tumorVolume= sio.loadmat(inputFilePath)
tree = tumorVolume[tumorVolumeNameInMat]
tree = (tree==1)

tumorPoints= sio.loadmat(inputFilePath)
treePoints = tumorPoints[tumorPointsNameInMat]

print("#######################################\n####### TISSUE CUTTER GENERATOR #######\n#######################################")
print("### Input file: "+inputFilePath)
print("##### Matrix name: "+tumorVolumeNameInMat)
print("### Output file: "+outputFilePath)

print("# Performing marching cubes algorithm on tumor volume")
verts, faces = measure.marching_cubes_classic(tree,)
tumorMesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        tumorMesh.vectors[i][j] = verts[f[j],:]
print("# done.")
print("# Save intermediate mesh")
tumorMesh.save("intermediate.stl")
print("# done.")

print("# MeshLab messsages:")
in_file = 'intermediate.stl'
out_file = 'intermediate_proc.stl'
filter_script_path = 'filter.mlx'
 # Add input mesh
command = "meshlabserver -i '"+os.getcwd()+"/" + in_file+"'"
# Add the filter script
command += " -s '"+os.getcwd()+"/" + filter_script_path+"'"
# Add the output filename and output flags
command += " -o '"+os.getcwd()+"/" + out_file + "' -om vn fn"
# Execute command

#print ("Going to execute: " + command)
output = subprocess.check_output(command, shell=True)
#last_line = output.splitlines()[-1]
#print(last_line)
#print("Done:")
#print (in_file + " > " + out_file)
#volume = ConvexHull(tree).volume
print("# done.")

print("# Convex hull generation for 2D projection:")
processed_tumour = mesh.Mesh.from_file('intermediate_proc.stl')
testArray = np.array([processed_tumour.vectors[:,0,0],processed_tumour.vectors[:,0,1]]).T
testArray.shape
tumorHull = ConvexHull(testArray,incremental=True)

if args.show_convex_hull:
    plt.figure(figsize=(10,10))
    plt.ylim([0,210])
    plt.xlim([0,210])
    plt.plot(processed_tumour.vectors[:,0,0],processed_tumour.vectors[:,0,1])
    plt.plot(testArray[tumorHull.vertices,0], testArray[tumorHull.vertices,1], 'r--', lw=2)
    plt.plot(testArray[tumorHull .vertices[0],0], testArray[tumorHull .vertices[0],1], 'ro')
    plt.show()


convexHullArray = np.array([[testArray[a,0], testArray[a,1]] for a in tumorHull.vertices])

convexHullPolygon =  Polygon([[testArray[a,0], testArray[a,1]] for a in tumorHull.vertices])
mbr_points = list(zip(*convexHullPolygon.minimum_rotated_rectangle.exterior.coords.xy))
mbr_lengths = [LineString((mbr_points[i], mbr_points[i+1])).length for i in range(len(mbr_points) - 1)]
minor_axis = min(mbr_lengths)
major_axis = max(mbr_lengths)



print("# done.")

print("# Build mould structure:")
convexHullExtrude = linear_extrude(slabHeight)(translate([0,0,0])(offset(r=5)(polygon([[testArray[a,0], testArray[a,1]] for a in tumorHull.vertices]))))
convexHullExtrude2 = linear_extrude(slabHeight)(translate([0,0,0])(offset(r=10)(polygon([[testArray[a,0], testArray[a,1]] for a in tumorHull.vertices]))))
tmpSlabs = []
numberOfSections = sizeInYDirection//sectionPlaneDistances if sectionDirection==1 else sizeInXDirection//sectionPlaneDistances
for i in range(numberOfSections):
    translateInX = i*sectionPlaneDistances if sectionDirection==0 else 0
    translateInY = i*sectionPlaneDistances if sectionDirection==1 else 0
    slabSize = [sizeInXDirection,sectionPlaneDistances-(gapBetweenSlabs/2),slabHeight] if sectionDirection==1 else [sectionPlaneDistances-(gapBetweenSlabs/2),sizeInYDirection,slabHeight]
    tmpSlabs += translate( [translateInX,translateInY,basePlateThickness])(
    cube(slabSize)
    )
d = intersection()(cube([sizeInXDirection,sizeInYDirection,basePlateThickness]),convexHullExtrude2)
d += intersection()(tmpSlabs,convexHullExtrude)
d -= translate( [0,0,basePlateThickness])(import_stl("intermediate_proc.stl"))

# Create vessel and contact point holes
for pointIdx,contactHole in enumerate(vesselPointIndices+contactPointIndices):
    centerOfMassCH = ndimage.measurements.center_of_mass(treePoints==contactHole)
    d = difference()(d,translate([centerOfMassCH[0],centerOfMassCH[1],0])(cylinder(h=sizeinZDirection,r=contactHoleSize)))
    d = difference()(d,translate([centerOfMassCH[0]+1.2*contactHoleSize,centerOfMassCH[1]+1.2*contactHoleSize,0])(mirror([1,0,0])(linear_extrude(1)(text(pointLabels[pointIdx], font = "Liberation Sans",size=5)))))


# Create tumor color part
centerOfMassTumor = ndimage.measurements.center_of_mass(treePoints==tumorPointIndex)
centerOfMassKidney = ndimage.measurements.center_of_mass(treePoints==kidneyPointIndex)

d_kidney = intersection()(translate([centerOfMassKidney[0],centerOfMassKidney[1],basePlateThickness-1])(cylinder(h=sizeinZDirection,r=contactHoleSize*1.5)),d)
d = difference()(d,translate([centerOfMassKidney[0],centerOfMassKidney[1],basePlateThickness-1])(cylinder(h=sizeinZDirection,r=contactHoleSize*1.5)))
d_tumour = intersection()(translate([centerOfMassTumor[0],centerOfMassTumor[1],basePlateThickness-1])(cylinder(h=sizeinZDirection,r=contactHoleSize*1.5)),d)
d = difference()(d,translate([centerOfMassTumor[0],centerOfMassTumor[1],basePlateThickness-1])(cylinder(h=sizeinZDirection,r=contactHoleSize*1.5)))


# Build a baseplate
# Determind minimum and maximum extend in alldirections
print(convexHullArray.shape)
minX = np.min(convexHullArray[:,0])
maxX = np.max(convexHullArray[:,0])
minY = np.min(convexHullArray[:,1])
maxY = np.max(convexHullArray[:,1])
print(minX,maxX,minY,maxY)
#print(sizeInXDirection/2-(maxX-minX)/2)
#d = translate([sizeInYDirection/2-(maxY-minY)/2,sizeInXDirection/2-(maxX-minX)/2,0])(d)
#d += translate([0,0,0])(cube([200,30,basePlateThickness]))


scad_render_to_file(d,'test.scad')
scad_render_to_file(d_tumour,'test_tumour.scad')
scad_render_to_file(d_kidney,'test_kidney.scad')
print("# done.")

#quit()
filesToConvert = ['test.scad','test_tumour.scad','test_kidney.scad']
targetFileNames = [outputFilePath,outputFilePath+"_tumour",outputFilePath+"_kidney"]
for fileIdx,fileToConvert in enumerate(filesToConvert):
    print("# Final conversion to .stl file:")
    # Add input mesh
    in_file = fileToConvert
    out_file = targetFileNames[fileIdx]+".stl"
    command = "openscad " + in_file+""
    # Add the output filename and output flags
    command += " -o " + out_file
    # Execute command

    print(command)
    print ("Going to execute: " + command)
    output = subprocess.check_output(command, shell=True)
    print("Done:")
    print (in_file + " > " + out_file)
    os.remove(fileToConvert)
#volume = ConvexHull(tree).volume
print("# Projected tumor size (max): \n## Major axis:"+str(major_axis)+" \n## Minor axis:"+str(minor_axis))
print("# done. Happy 3D printing!")
# tidy up
os.remove("intermediate.stl")
os.remove("intermediate_proc.stl")
