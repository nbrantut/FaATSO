function [Nx,Ny,Nz,h,data] = import_data(filename)
%[Nx,Ny,Nz,h,data] = IMPORT_DATA(filename)
%
%this function imports an array "dat" of size Nx-by-Ny-by-Nz, and a grid
%spacing parameter "h", from a binary file "filename" which is an output
%from fmm.
%
%input:
%   filename:   name of binary file in which the data will be written.
%
%output:
%   Nx:         number of points in x direction
%   Ny:         number of points in y direction
%   Nz:         number of points in z direction
%   h:          grid spacing
%   data:       an Nx-Ny-Nz array with the data.

fid = fopen(filename);
N = fread(fid,3,'int');
h = fread(fid,1,'double');
T = fread(fid,'double');
fclose(fid);

Nx = N(1);
Ny = N(2);
Nz = N(3);


data = reshape(T,[Nx Ny Nz]);
