function [Vh_prior, E_prior, events_prior] = priormodel_2d_xy(V0, Vh_true, events_true, noise_p)
%[Vh_prior, E_prior, events_prior] = PRIORMODEL_2D_XY(V0, Vh_true, events_true, noise_p)
%
%This function builds a homogeneous, isotropic prior velocity model as well
%as "noisy" prior event locations based on the true ones.
%
%input:
%    V0:          scalar value of the homogeneous prior velocity
%    Vh_true:     the true Vh model (size Nx Ny Nz)
%    events_true: the true event locations of the synthetic test
%    noise_p:     the noise added on the event positions to get the prior
%
%output:
%    Vh_prior:    the prior velocity model
%    E_prior:     the prior anisotropy model (0 everywhere)
%    events_prior:the prior noisy event locations

%initialise prior Vh
Vh_prior = V0 + 0*Vh_true;
%keep the "mask" of zero values outside the disk of radius 2.5:
Vh_prior(Vh_true==0) = 0;

%initialise prior E
E_prior = 0*Vh_true;

%get number of events
Nev = size(events_true, 1);

% add noise to event positions
events_prior = events_true;
events_prior(:,1) = events_prior(:,1) + (rand(Nev,1)*2 - 1)*noise_p;
events_prior(:,2) = events_prior(:,2) + (rand(Nev,1)*2 - 1)*noise_p;
events_prior(:,4) = events_prior(:,4) + (rand(Nev,1)*2 - 1)*noise_p;


