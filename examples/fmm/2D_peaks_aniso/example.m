%2D_peaks_aniso example
%
%This is a basic test of the fast marching method in 2d in a heterogeneous,
%isotropic medium. This Matlab program runs the program fmm using a
%square grid of 201x201 nodes, with spacing h=0.1 mm and a velocity and 
% anisotropy given by Matlab's "peaks" function, with a source located at 
%random. 
%
%The arrival times are displayed as contour, and a small animation of the
%wavefront propagation is produced.
%
%HOW TO RUN THIS EXAMPLE:
%
% 1. make sure that the wrapper function run_fmm is in your Matlab path.
% 2. maybe modify the path to the fmm binary if needed
% 3. run the script.

%set path to fmm binary
%NOTE: change this if you decided to build fmm with a different name or
%path!
fmmpth = '../../../bin/./fmm';

%spacing
h =.1;

%number of points in each direction
Nx = 201;
Ny = 1;
Nz = 201;

%initialise V and epsilon

Vh = reshape(3 + peaks(Nx)/3, Nx, 1, Nz);
E =  reshape(peaks(Nx)/20, Nx, 1, Nz);

%coordinates of the source in real dimensions
source = rand(1,3).*[Nx Ny Nz].*h;

%say yes we want a source box
box = 1;

% compute arrival times using fmm
[T,x,y,z] = run_fmm(source, Vh, E, h, box, fmmpth);

%plot result in contourplot
show_movie_2d_xz(x,z,Vh,E,T, 0.01)