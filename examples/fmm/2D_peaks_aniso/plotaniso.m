function plotaniso(x,z,E,step)
%PLOTANISO(x,z,E,step)
%
%This function plots a thin grey line at regularly spaced nodes withi a 2D
%x-z grid. The orientation of the line gives the direction of fast P-wave,
%and the length of the line gives the relative difference between fast and
%slow velocity.
%
%input:
%    x:    x coordinates (1d array)
%    z:    z coordinates (1d array)
%    E:    anisotorpy (Nx,1,Nz) array
% step:    step between lines (if this is too small, this function will be
% very slow and display too many lines... use something around 10-20.)

%retrieve grid dimensions
h=x(2)-x(1);
Nx=length(x);
Nz=length(z);

for ix=1:step:Nx
    for iz=1:step:Nz
        if E(ix,1,iz)>0
            plot(x(ix)+[0 0], z(iz)+[-h/2 h/2]*step*(E(ix,1,iz)),'color', [0.3 0.3 0.3]);
        else
            plot(x(ix)+[-h/2 h/2]*step*(-E(ix,1,iz)./(1+E(ix,1,iz))), z(iz)+[0 0],'color', [0.3 0.3 0.3]);
        end
    end
end