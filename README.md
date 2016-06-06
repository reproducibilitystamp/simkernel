# simkernel
This code reproduces the results in the following paper,

Complex Transfinite Barycentric Mappings with Similarity Kernels
<br>Comput. Graph. Forum. (Proc. Sym. Geom. Proc.), 35(5), 2016
<br>By    Renjie Chen and Craig Gotsman

# Requirements
- MS Windows (Windows 7/8/10)
- MATLAB (>2014a)

# How to run the code
Execute the script fig4.m in MATLAB, a figure will be generated, which will be a 
reproduction of figure 4 in the paper.
Similarly, to reproduce figure 5 in the paper, simply run fig5.m.

# What the code does
The script will load the source contour and cubic spline parameters for generating the target
contour from file (fig4_data.mat and fig5_data.mat). The source contour will be subsampled 
and triangulated using the triangle https://www.cs.cmu.edu/~quake/triangle.html software. The
transfinite mean value coordiantes mapping will be applied to each vertex of the triangle
mesh (source mesh) to the target contour, then the iterative inverse mapping algorithm in 
Section 7 of the paper will be applied to the vertices and recover their positions in the 
source mesh.
