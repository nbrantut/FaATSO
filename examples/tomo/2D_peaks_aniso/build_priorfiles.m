function build_priorfiles(Vh_prior, E_prior, h, events_prior, folder)
%build_priorfiles(Vh_prior, E_prior, events_prior, folder)
%
%This function generate files for prior velocity, anisotropy and event
%location models.
%
%input:
%    Vh_prior:     the prior velocity model
%    E_prior:      the prior anisotropy model
%    h:            grid spacing
%    events_prior: the prior event locations
%    folder:       folder where to store the files

[Nx,Ny,Nz] = size(Vh_prior);

%write events file
dlmwrite([folder  'events_prior.txt'],events_prior,' ');

% export V and E models
export_data([folder 'V_prior.bin'], Nx, Ny, Nz, h, log(Vh_prior));
export_data([folder 'E_prior.bin'], Nx, Ny, Nz, h, E_prior);