function M = m2fields(m, Nx, Ny, Nz)
%M = M2FIELDS(m, Nx, Ny, Nz)
%
%This function extracts individual fields (velocity, anisotropy, AE
%locations) from the assembled vector m of model parameters as used by the
%tomography program.
%
%input:
%    m:   the vector of model parameters
%    Nx, Ny, Nz: the grid dimension
%
%output:
%    M:  structure with fields:
%         .V0:  an Nx*Ny*Nz array with the wave velocity
%         .E0:  an Nx*Ny*Nz array witht the anisotropy parameter
%         .N:   a vector of North coordinates of the sources
%         .E:   a vector of East coordinates of the sources
%         .D:   a vector of Depth coordinates of the sources
%         .t0:   a vector of origin times of the sources

Ntot = Nx*Ny*Nz;

M.V0 = exp(reshape(m(1:Ntot), [Nx Ny Nz]));
M.E0 = reshape(m(Ntot+1:2*Ntot), [Nx Ny Nz]);

M.N  = m(2*Ntot+1:4:end);
M.E  = m(2*Ntot+2:4:end);
M.D  = m(2*Ntot+3:4:end);
M.t0 = m(2*Ntot+4:4:end);
