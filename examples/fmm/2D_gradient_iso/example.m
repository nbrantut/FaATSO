%2D_gradient_iso
%
%This routine computes a benchmarking test of the second order fmm, based
%on a known analytical solution for rays propagating in a medium with
%velocity linearly increasing with depth.
%
%This test basically replicates the figure 10c of Rawlinson and Sambridge,
%GJI 2004, using similar parameters and the same methods.
%
% the source is at the surface, at position X=0.
% the medium has a velocity of 4 mm/us at the surface, and it increases by
% 0.1 mm/us every mm (gradient of 0.1 /us).
%
% four grid spacing are tested: 0.125, 0.25, 0.5 and 1 mm spacing.
%
%HOW TO RUN THIS EXAMPLE:
%
% 1. make sure that the wrapper function run_fmm is in your Matlab path.
% 2. make sure that the executable fmm is in the PATH (check: getenv('PATH'))
% 3. run the script.

%% generate grid

%spacing
tab_h =[1  0.5 0.25 0.125];

%number of points in each direction
tab_Nx = [41 81 161 321];
tab_Ny = [101 201 401 801];
Nz = 1;

%% compute solution for all grids:

%preallocate array of strcut. to store results:
sol = repmat(struct('y',{}, 'T', {}, 'T_true', {}),4,1);

for it=1:4
    
    h = tab_h(it);
    Nx = tab_Nx(it);
    Ny = tab_Ny(it);
    
    %'vertical' coord:
    x = h*(0:Nx-1);

    %initialise V and epsilon
    Vh = 4 + zeros(Nx,Ny,Nz);
    E = 0.0 + zeros(Nx,Ny,Nz);

    for i=2:Nx
        Vh(i,:,:) = 4 + 0.1*x(i);
    end

    % source point

    %coordinates of the source in real dimensions
    source = [0,0,0];
    
    box = 0;

    % compute arrival times
    [T,x,y,z] = run_fmm(source, Vh, E, h, box);
    
    % extract only what will be useful for our comparison
    sol(it).y = y;
    sol(it).T = T(1,:,1);
    sol(it).T_true = Tsurf_Vlineardepth(4,0.1,y);
    
end


%% make the plots

%first plot: trace wavefronts using the most accurate version, for
%display:
figure;
contour(y,x, T, 10, 'k');
axis equal;
set(gca, 'ydir', 'reverse',...
    'xlim', [y(1)  y(end)],...
    'ylim', [x(1)  x(end)]);
xlabel('{\itx} (mm)');
ylabel('{\ity} (mm)');


%comparison between grid sizes, plot similar to fig 10c of Rawlinson and Sambridge, GJI 2004
figure;
plot(sol(4).y, 1e3*abs(sol(4).T - sol(4).T_true));
hold all
plot(sol(3).y, 1e3*abs(sol(3).T - sol(3).T_true));
plot(sol(2).y, 1e3*abs(sol(2).T - sol(2).T_true));
plot(sol(1).y, 1e3*abs(sol(1).T - sol(1).T_true));

set(gca,'xlim',[0 100],...
    'ylim',[0 80]);

xlabel('{\itx} (mm)');
ylabel('|{\itt} - {\itt}_{true}| (ns)');

legend(['h='  num2str(tab_h(4)) ' mm'],...
    ['h='  num2str(tab_h(3)) ' mm'],...
    ['h='  num2str(tab_h(2)) ' mm'],...
    ['h='  num2str(tab_h(1)) ' mm']);
