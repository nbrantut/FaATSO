%2D_homogeneous_iso example
%
%This is a basic test of the fast marching method in 2d in a homogeneous,
%isotropic medium. This Matlab program runs the program fmm using a
%square grid of 201x201 nodes, with spacing h=0.1 mm and a velocity of 1
%mm/us, with a source located at one corner. Then the arrival time are
%imported and displayed as contours, together with the analytical solution
%for comparison.
%
%HOW TO RUN THIS EXAMPLE:
%
% 1. make sure that the wrapper function run_fmm is in your Matlab path.
% 2. make sure that the executable fmm is in the PATH (check: getenv('PATH'))
% 3. run the script.

%spacing
h =.1;

%number of points in each direction
Nx = 201;
Ny = 201;
Nz = 1;

%initialise V and epsilon
v0 = 1;
e0 = 0;

Vh = v0* ones(Nx);
E =  e0 *ones(Nx,Ny,Nz);

%coordinates of the source in real dimensions
source = [0 0 0];

%say yes we want a source box
box = 1;

% compute arrival times using fmm
[T,x,y,z] = run_fmm(source, Vh, E, h, box);

% compute exact arrival times using analytical formula
xx = meshgrid(x,y);
yy = meshgrid(y,x)';
Texact = sqrt((xx-source(1)).^2 + (yy-source(2)).^2)./v0;

%plot result in contourplot
figure;

contour(x,y,T',10,'r');
hold on
contour(x,y,Texact,10,'k--');
axis equal;
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
xlabel('{\itx} (mm)');
ylabel('{\ity} (mm)');
title('{\itV}_h = 1 mm/\mus, {\itE} = 0');