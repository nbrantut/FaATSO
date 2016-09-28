function export_data(filename, Nx,Ny,Nz,h,data)
%export_data(filename, dat,Nx,Ny,Nz,h)
%
%this function takes an array "dat" of size Nx-by-Ny-by-Nz, and a grid
%spacing parameter "h", and writes it in a binary file "filename" for use
%as input in the fmm code
%
%input:
%   filename:   name of binary file in which the data will be written.
%   Nx:         number of points in x direction
%   Ny:         number of points in y direction
%   Nz:         number of points in z direction
%   h:          grid spacing
%   data:       an Nx-Ny-Nz array with the data.
%
%output: NONE;

fid = fopen(filename,'w');

fwrite(fid,Nx,'int');
fwrite(fid,Ny,'int');
fwrite(fid,Nz,'int');
fwrite(fid,h,'double');

VV = data(:);

fwrite(fid, VV, 'double');

fclose(fid);
