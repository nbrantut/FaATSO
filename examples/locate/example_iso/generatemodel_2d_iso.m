function [Vh, E, h, sensors, events] = generatemodel_2d_iso(Nev)
%[Vh, E, h, sensors, events] = GENERATEMODEL_2D_ISO(Nev)
%
%This function constructs a 2D model of 28x28 nodes with an isotropic,
%heterogeneous wave velocity (given by matlab's peaks function), and
%generates coordinates for 8 sensors positioned regularly in the grid on
%the perimeter of a disk of radius 2.5. OUtisde this disk, the velocity is
%set to 0. The function then generates a set of Nev random source positions
%inside the disk.
%
%input:
%    Nev:    number of events
%
%output:
%    Vh:      the velocity structure, given by a scaled "peaks" function,
%             and set to 0 outside a disk of radius 2.5.
%    E:       the anisotropy strucutre (0 everywhere).
%    h:       the grid spacing (=0.2);
%    sensors: array of sensor positions
%    events:  array of event positions


%% generate grid with Vh and E fields

%spacing
h=2;

%number of points in each direction
Nx=28;
Ny=28;
Nz=1;

%Coordinates
x = h*(0:Nx-1);
y = h*(0:Ny-1);

%initialise V and epsilon
Vh = reshape(3 + peaks(Nx)/10, [Nx Ny Nz]);

XX = meshgrid(x,y);
YY = meshgrid(y,x)';
x_center = x(end)/2;
y_center = y(end)/2;
RR = sqrt((XX-x_center).^2 + (YY-y_center).^2);
ind_out = find(RR>28);
[Iout,Jout] = ind2sub([Nx Ny],ind_out);

for k=1:length(Iout)
    Vh(Iout(k),Jout(k),:) = 0.0;
end

E = zeros(Nx,Ny,Nz);


%% generate sensors and events

%radius positions of sensors:
r_sens = 25;
%azimuths of sensors:
q_sens = pi/4 * (0:7);
%build sensor locations matrix:
sensors = repmat([x_center  y_center  0  0],8,1) + ...
    r_sens.*[cos(q_sens')  sin(q_sens')  0*q_sens'  0*q_sens'];

%pick random radii for events between 0 and 2.5
r_evt = rand(Nev, 1)*25;
%pick random azimuths for eveents between 0 and 2pi
q_evt = rand(Nev, 1)*2*pi;
%build event locations matrix
events = [r_evt.*cos(q_evt)+25+h  r_evt.*sin(q_evt)+25+h zeros(Nev,1) rand(Nev,1)*4];
