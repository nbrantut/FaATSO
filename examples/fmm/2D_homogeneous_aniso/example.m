%2D_homogeneous_aniso example
%
%This is a basic test of the fast marching method in 2d in a homogeneous,
%anisotropic medium. This Matlab program runs the program fmm using a
%square grid of 201x201 nodes, with spacing h=0.1 mm and a velocity of 1
%mm/us for an aniostropy coefficient E=0.25, with a source located at one 
%corner. Then the arrival time are imported and displayed as contours, 
%together with the analytical solution for comparison.
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
Vh = ones(Nx,Ny,Nz);
E = 0.25 + zeros(Nx,Ny,Nz);

%coordinates of the source in real dimensions
source = [0 0 0];

%say yes we want a source box
box = 1;

% compute arrival times
[T,x,y,z] = run_fmm(source, Vh, E, h, box, fmmpth);

% compute exact arrival times using analytical formula
Texact = analytical_2d_hom_aniso(source, x, z, 1, 0.25);

%plot result in contourplot
contour(x,z,squeeze(T)',10,'r');
hold on;
contour(x,z,Texact',10,'k--');
axis equal;
xlim([x(1) x(end)]);
ylim([z(1) z(end)]);
xlabel('{\itx} (mm)');
ylabel('{\itz} (mm)');
title('{\itV}_h = 1 mm/\mus, {\itE} = 0.25');