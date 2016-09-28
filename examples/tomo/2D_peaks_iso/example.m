%2D_peaks_iso example
%
%This is a basic test of the tomographic method, using a 2d isotropic
%heterogeneous velocity field, 8 transducers and a given number of events
%as passive sources.
%
%The true model is computed using a 28x28 grid with grid spacing 0.2. Only 
%a disk of radius has a non zero velocity (this is to test the method when
%using a "mask" of V=0 in a given domain). A set
%number of events are generated at random within the disk, and 8 sensors 
%are positioned around it.
%
%The true arrival times between events and sensors, and between pairs of
%sensors, are computed using the program "synthetic". A set level of
%artifical noise is added to the true arrival times.
%
%The prior model is initially homogeneous, isotropic, and the prior event
%positions are the true ones plus a set level of random "noise".
%
%The inversion problem is solved using the "tomo" program, and the results
%are displayed: one plot shows the inverted mean V model, the other shows
%the true model, and the true and inverted event positions are plotted as
%well. The residual vs. iteration number is also shown.
%
%HOW TO RUN THIS EXAMPLE:
%
% 1. mayeb modify the path to the "synthetic" binary if needed.
% 2. make sure that the wrapper function run_synthetic is in your Matlab path.
% 3. maybe modify the path to the tomo binary if needed
% 4. run the script.
%
%You may want to modify the parameters.txt file to explore the effect of
%changing inversion parameters etc.


% number of events to try:
Nev = 40;

% noise level on arrival times
noise_t = 0.05;

% prior homogeneous V value
V0 = 3;

% noise on prior locations
noise_p = 0.2;

% name of folder where to store the files generated:
folder = 'files/';

% path to the executable to run the synthetic data
synth_path = '../synthetic/./synthetic';
% path to the run_synthetic wrapper function:
addpath('../synthetic');

% path to the tomo binary
tomo_path = '../../../bin/./tomo';

%first check that the folder exists, otherwise create it:
if (exist(folder,'dir')~=7)
    mkdir(folder);
end
    
% compute true model with the set number of events:
[Vh, E, h, sensors, events] = truemodel_2d_xy(Nev);

% compute synthetic data with noise:
run_synthetic(events, sensors, Vh, E, h, noise_t, synth_path, folder)

% compute prior model
[Vh_prior, E_prior, events_prior] = priormodel_2d_xy(V0, Vh, events, noise_p);

% make files from prior model
build_priorfiles(Vh_prior, E_prior, h, events_prior, folder);

% run tomo:
cmd = [tomo_path  ' parameters.txt'];
system(cmd);

%% show results in a plot:
display_results_xy;

%% show posterior movie:
posterior_movie_xy('postmovie.gif', ...
    20, ...
    [2.6  3.4],...
    200, ...
    [folder 'mpost.txt'], ...
    [folder 'CMpost.arma.bin'],...
    Nx, Ny, Nz, h);