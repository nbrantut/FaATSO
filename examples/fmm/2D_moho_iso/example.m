%2D_moho_iso example
%
%This is a basic test of the fast marching method in 2d in a two-layer,
%isotropic medium. This Matlab program runs the program fmm using a
%square grid of 201x801 nodes, with spacing h=0.1 mm and two layers with 
%a velocity of 1mm/us at the top and 4 at the bottom, with a source 
%located at one corner. Then the arrival time are imported and displayed
%as contours, together with the surface hodochrones.
%
%Note that the velocity contrast is indeed a very strong one (much
%stronger than the typical change at the moho...), but this is for
%illustrative purposes only.
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
Ny = 801;
Nz = 1;

%initialise V and epsilon
Vh = ones(Nx,Ny,Nz);
E =  zeros(Nx,Ny,Nz);

%set fast layer at the bottom
i_moho = 121;
Vh(i_moho:end,:,:) = 4;


%coordinates of the source in real dimensions
source = [0 0 0];

%say yes we want a source box
box = 1;

% compute arrival times
[T,x,y,z] = run_fmm(source, Vh, E, h, box);

%plot result in contourplot
figure;
contour(y,x,T,10,'k');
hold on
plot(y,y*0 + x(i_moho),'k:')
axis equal;
set(gca, 'xlim',[y(1) y(end)],...
    'ylim',[x(1) x(end)],...
    'ydir','reverse');
xlabel('{\itx} (mm)');
ylabel('{\ity} (mm)');
title('{\itV}_h = 1 to 4 mm/\mus, {\itE} = 0');

%plot hodochrones
figure;
plot(y, T(1,:,1));
hold on;
plot(y,y,'k:');

d = (i_moho-1)*h;
plot(y,2*d*sqrt(1-(1/4)^2) + y/4, 'k:');
xlabel('position from source (mm)');
ylabel('arrival time (\mus)');
set(gca, 'xlim', [y(1)  y(end)],...
    'ylim', [0, max(T(1,:,1))])