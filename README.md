## 3D printed tumor moulds for guided tissue sampling
This software should be cited with the DOI 10.5281/zenodo.3066305.

## Requirements
* MATLAB
* CERR (https://github.com/cerr/CERR)
* Python (with modules as listed in requirements.txt)
* Meshlab (with meshlabserver)
* OpenSCAD

## Examples

### Rotation and habitat map generation
Images are processed in MATLAB to generate habitat maps and correctly oriented volumes for printing. To test the code, run the following commands in a MATLAB prompt:
```
>> rotateTumour /path/to/CERR
```
This script creates two .mat files: 
1. transf_patient.mat (containing the transformation matrices required to rotate the image to the desired orientation)
2. aligned_tumour.mat (containing a 3D matrix with the volume to print, already rotated)
In addition, the script creates multiple 2D and 3D visualizations of the tumour.

```
>> makeHabitats /path/to/CERR
```
This script creates one .png image per mould slice, containing the habitat maps. In addition, it creates a series of histograms comparing multiple imaging parameters for each habitat.

### Mould generation
`python generate_mould.py example_data/aligned_tumour.mat`
will result in three .stl files:
1. Patient.stl (containing the mould with holes for projected centroids of tumour and kidney)
2. Patient_tumour.stl (containing the fill-in for the tumour centroid hole)
3. Patient_kidney.stl (containing the fill-in for the kidney centroid hole) cutter
