function T_anal = analytical_2d_hom_aniso(source, x, z, Vh, E)
%T_anal = ANALYTICAL_2D_HOM_ANISO(source, x, z, Vh, E)
%
%This function computes the analytical solution for arrival times in a 2D,
% homogeneous anistropic medium.
%
%input:
%    source:   source position (x,0,z)
%    x:        vector of x coord.
%    z:        vector of z coord.
%    Vh:       horizontal velocity
%    E:        anisotropy parameter
%
%ouptut:
%    T_anal:   array of arrival times from source.

%initialise matrix:
T_anal = zeros(length(x),length(z));

%source coord:
x0 = source(1);
z0 = source(3);


for ix=1:length(x)
    for iz=1:length(z)
        %group angle:
        phi = atan((x(ix)-x0)/(z(iz)-z0));
        if isfinite(phi)
            %group velocity:
            V_phi = group_vel(phi, Vh, E);
            %arrival time:
            T_anal(ix,iz) = sqrt( (x(ix)-x0)^2 + (z(iz)-z0)^2)./V_phi;
        end
    end
end