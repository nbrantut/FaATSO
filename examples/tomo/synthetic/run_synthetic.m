function run_synthetic(events, sensors, Vh, E, h, noise, folder)
%RUN_SYNTHETIC(events, sensors, Vh, E, h, noise, synth_path, folder)
%
%This functions build files for synthetic dataset, and stores then into a
%specified folder, which must exist prior to calling the function (otherwise it will throw an error).
%
%NOTE: this function will return an error if the program "synthetic" is not in
%the PATH environment variable used by Matlab. To fix this, use the routine
%"put_synthetic_on_path" or directly add the path to fmm executable using
%setenv.
%
%input:
%    events:    array of event positions
%    sensors:   array of sensor positions
%    Vh:        velocity model
%    E:         anisotropy model
%    h:         grid spacing
%    noise:     noise level to add to the arrival times
%    folder:    folder where to store the output files

[Nx,Ny,Nz] = size(Vh);

if any(size(Vh)~=size(E))
    disp('Input velocity and anisotropy grids are not compatible.');
    return;
end

sensfile = [folder  'sensors_synth.txt'];
evtfile = [folder  'events_synth.txt'];

dlmwrite(sensfile,sensors,' ');
dlmwrite(evtfile,events,' ');

Vfile = [folder  'V_synth.bin'];
Efile = [folder  'E_synth.bin'];

export_data(Vfile, Nx, Ny, Nz, h, log(Vh));
export_data(Efile, Nx, Ny, Nz, h, E);

%% run synthetic test
tshotfile = [folder  't_shots_synth.txt'];
tevtfile = [folder  't_events_synth.txt'];

cmd = ['synthetic '  num2str(noise) ' ' Vfile ' ' Efile ' ' sensfile ' ' evtfile ' ' tshotfile ' ' tevtfile];

system(cmd);