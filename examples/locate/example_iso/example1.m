%2D_peaks_iso example
%
%This is a basic test of the location method, using a 2d isotropic
%heterogeneous velocity field, 8 transducers and a given number of events
%as passive sources.
%
%The true model is computed using a 28x28 grid with grid spacing 2. Only 
%a disk of radius 25 has a non zero velocity (this is to test the method when
%using a "mask" of V=0 in a given domain). A set
%number of events are generated at random within the disk, and 8 sensors 
%are positioned around it.
%
%The true arrival times between events and sensors, and between pairs of
%sensors, are computed using the program "synthetic". A set level of
%artifical noise is added to the true arrival times.
%



% number of events to try:
Nev = 40;

% noise level on arrival times
noise_t = 0.05;

% prior homogeneous V value
V0 = 3;

% name of folder where to store the files generated:
folder = 'files/';

%first check that the folder exists, otherwise create it:
if (exist(folder,'dir')~=7)
    
    mkdir(folder);
end
    
% compute true model with the set number of events:
[Vh, E, h, sensors, events] = generatemodel_2d_iso(Nev);

% save in files 
[Nx,Ny,Nz] = size(Vh);

export_data([folder 'vh.bin'],Nx, Ny, Nz, h, log(Vh));
export_data([folder 'E.bin'],Nx, Ny, Nz, h, E);
dlmwrite([folder 'sensors.txt'],sensors,' ');
dlmwrite([folder 'events.txt'],events,' ');

% compute synthetic data with noise:
!synth_loc params_synth.txt

% get locations from least L1 norm grid search
!loc params.txt

%% plot outcome
%import posterior event locations
events_out = dlmread('files/events_post.txt');

%plot results
figure;

pcolor(h*(0:Nx-1), h*(0:Ny-1),Vh');
axis equal;
xlim([0 54]);
ylim([0 54]);
caxis([2 4])

hold on
plot(sensors(:,1), sensors(:,2),'v','markersize',8);
plot(events(:,1), events(:,2),'k+');
plot(events_out(:,1), events_out(:,2),'ro');
