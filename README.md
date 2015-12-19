# Dense Histogram Computation on Polygonal Regions
Compute histograms densely on polygonal regions

Code for the following paper:

> Yuting Zhang, Yueming Wang, Gang Pan, Zhaohui Wu, “Efficient Computation of Histograms on Densely Overlapped Polygonal Regions”, *Neurocomputing*, vol. 118, pp 141 – 149, 22 October 2013. 

All the methods addressed in this paper are included in this code. 

Dependencies
===

No particular dependency for the library

OpenCV Library (for demo)

Compilation
===
Use `gcc`:

	make

Use `icpc`:

	make -f Makefile.icc
	
Use `cygwin`:

	make -f Makefile.cygwin
	
Manual
===

CMD: `polyint-static-demo`

ARGS:
 
- `-m`  method
	- `-2` - integral image
	- `-1` - naive
	- `0` - fast auto
	- `1` - fast edge
	- `2` - fast naive
- `-f`  input image path
- `-b`  binnum base, binnum = N^2
- `-p`  polygon file (Can use multiple times for multiple polygon)
	example: polyint/polygon.txt (specifying vertices)
- `-z`  polygon zoom factor
- `-t`  function type
	- `1` - linear
	- `2` - Chi
	- `3` - Entropy
