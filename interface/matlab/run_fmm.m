function [T,x,y,z] = run_fmm(source, Vh, E, h, box)
%[T,x,y,z] = RUN_FMM(source, Vh, E, h, box, fmm_path)
%
%This function is a wrapper to run the fmm program based on source location
%and velocity structures defined in matlab.
%
%NOTE: this function will return an error if the program "fmm" is not in
%the PATH environment variable used by Matlab. To fix this, use the routine
%"put_faatso_on_path" or directly add the path to fmm executable using
%setenv.
%
%input:
%    source:   a triplet of (x,y,z) coordinates where the source will be.
%    Vh:       an Nx-by-Ny-by-Nz array containing values of the horizonal
%              velocity at each node.
%    E:        an Nx-by-Ny-by-Nz array containing values of the anisotropy
%              coefficient at each node.
%    h:        the grid spacing.
%    box:      0/1 depending if a source box is to be used.
%
%output:
%    T:    an Nx-by-Ny-by-Nz array containing values of the arrival time
%          at each node.
%    x:    a mesh with x coordinates of the nodes
%    y:    a mesh with y coordinates of the nodes
%    z:    a mesh with z coordinates of the nodes

%% generate grid

[Nx,Ny,Nz] = size(Vh);

if any(size(Vh)~=size(E))
    disp('Input velocity and anisotropy grids are not compatible.');
    T=[];
    x=[];
    y=[];
    z=[];
    return;
end

%Coordinates 
x = h*(0:Nx-1);
y = h*(0:Ny-1);
z = h*(0:Nz-1);


%% source point
x0 = source(1);
y0 = source(2);
z0 = source(3);

%find nearest grid point and extract its index
[~,i] = min(abs(x0 - x));
[~,j] = min(abs(y0 - y));
[~,k] = min(abs(z0 - z));

ind0 = (i-1) + (j-1)*Nx + (k-1)*Nx*Ny; 


%% export V and E models
safename = randstr(25);
Vfile = [safename '_V.bin'];
Efile = [safename '_E.bin'];

export_data(Vfile, Nx, Ny, Nz, h, Vh);
export_data(Efile, Nx, Ny, Nz, h, E);

%% run the fmm algorithm
Tfile = [safename '_T.bin'];
cmd = ['fmm ' Vfile ' ' Efile ' ' num2str(ind0) ' ' num2str(box) ' ' Tfile];
system(cmd);


%% import the results
[~,~,~,~,T] = import_data(Tfile);

%% clean the directory
cmd = ['rm ' Vfile ' ' Efile ' ' Tfile];
system(cmd);
