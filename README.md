## 3D printed tumor moulds for guided tissue sampling
## Requirements
* MATLAB
* CERR (https://github.com/cerr/CERR)
* Python (with modules as listed in requirements.txt)
* Meshlab (with meshlabserver)
* OpenSCAD

## Examples

### Rotation and habitat map generation
Images are processed in MATLAB to generate habitat maps and correctly oriented volumes for printing. To test the code, run the following commands in a MATLAB prompt:
    >> rotateTumour /path/to/CERR
    >> makeHabitats /path/to/CERR

### Mould generation
`python generate_mould.py example_data/aligned_tumour.mat`
will result in three .stl files:
1. Patient.stl (containing the mould with holes for projected centroids of tumour and kidney)
2. Patient_tumour.stl (containing the fill-in for the tumour centroid hole)
3. Patient_kidney.stl (containing the fill-in for the kidney centroid hole) cutter
