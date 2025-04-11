# ConcaveCoverSolver
Introduction to this project
	This program is used to automate the extraction of polygon points for any arbitrary single graph, and the convex split of concave polygons, using automated scripts to generate data file Agen.f90 that can be used as input to a thinnest coverage problem (TCP) solving toolbox. In this project, polygon vertex extraction, convex cutting and data file generation were carried out for square, circle, Minkowski island fractal and a hand-painted graph, and the TCP problem solving was carried out with the toolbox.

Folder description
	modA_auto_generation: MATLAB code that automatically generates the input data file Agen.f90 for TCP-solving toolbox.
	modified_files_in_bglcovering-main: Fortran code and bash code used to modify files with the same name in the TCP-solving toolbox.
	TCP_results: TCP solution results for four graphs with 1000 trials.

Running environment
	modA_auto_generation：Windows 11, MATLAB R2019a
	modified_files_in_bglcovering-main: Ubuntu GNU/Linux (version 20.04.3 LTS), GFortran compiler of GCC (version 9.4.0)

Dependent file
（1）TCP-solving toolbox
	 website: https://github.com/johngardenghi/bglcovering.
	 reference: [1] E. Birgin, J. L. Gardenghi, and A. Laurain, “Asymptotic bounds on the optimal radius when covering a set with minimum radius identical balls ∗,” 2022. Accessed: Mar. 20, 2025.
			   [2] E. G. Birgin, J. L. Gardenghi, and A. Laurain, “Bounds on the Optimal Radius When Covering  a Set with Minimum Radius Identical Disks,” Mathematics of OR, vol. 49, no. 3, pp. 1855–1889, Aug. 2024.
	  
	 For random initial guesses rather than the lattice-based initial guesses in [1], also see 
	 website: https://www.ime.usp.br/~egbirgin/sources/coveringsecond/
	 reference: [3] E. G. Birgin, A. Laurain, R. Massambone, and A. G. Santana, “A Shape-Newton approach to the problem of covering with identical balls,” SIAM J. Sci. Comput., vol. 44, no. 2, pp. A798–A824, Apr. 2022.
  (2)   The third party codes
	 Algencan 4.0.0
	 website: https://www.ime.usp.br/~egbirgin/sources/bmcomper/

	 GEOMETRY 
	 website: https://people.sc.fsu.edu/~jburkardt/f_src/geometry/geometry.html
 
	 GEOMPACK
	 website: https://people.math.sc.edu/Burkardt/f_src/geompack2/geompack2.html
       
         BLAS
         type in terminal:
	 	sudo apt update
	  	sudo apt install libblas-dev
	 HSL
         website: https://licences.stfc.ac.uk/products
         In this project, LibHSL - 2023.11.7 is adopted.

Running steps:
（1）Install the aforementioned toolbox and dependent libraries.
（2）In .\modA_auto_generation\polygon_extract.m,
		Select input of graphic data: img = imread('Minkowski_island_fractal.png');
		Select the boundary: contour = boundaries{2};
		Set the number of smoothing point: x= smooth(x,3)；y = smooth(y,3);
		Set the threshold for detecting initial high curvature points: threshold = k3(3)*1.0;
		Set the maximum pixel point error: max_error =5.0;
	Then, run the program to generate polygon data.
（3）In .\modA_auto_generation\modA_autogen.m
		Select polygon data: load MIF_polygon_parameters.mat;
	Then, run the program to generate the data file Agen.f90 after polygon split.
（4）Replace the files with the same names in the TCP-solving toolbox's running path with those in the ./modified_files_in_bglcovering-main path.
（5）Replace the Agen.f90 file in the running path folder with the corresponding Agen.f90 file in (3).
（6）Enter the command in the terminal under the running path：./run-tests.sh Agen random
（7）In the output file named "output-....dat", r represents reference radius and pol area represents reference area.
