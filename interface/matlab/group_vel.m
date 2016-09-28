function Vg = group_vel(phi, Vh, E)
%Vg = GROUP_VEL(phi, Vh, E)
%
%This function computes the group velocity based on the input group angle
%(i.e., the angle formed by the ray direction and the vertical axis), the
%horizontal velocity and the anisotropy parameter, using the simple
%elliptical anisotropy model for phase velocity:
%
%  V_phase (theta)  = Vh * (1 + E*cos(theta)^2 )
%
% where theta is the phase angle.
%
%input:
%    phi:    group angle
%    Vh:     horizontal velocity
%    E:      anisotropy parameter
%
%output:
%    Vg:     group velocity

%determing phase angle based on input gropu angle:
theta = fzero(@(th) group_angle(th,E) - phi, phi);

%group velocity computation:
Vg = Vh*(1 + E*cos(theta).^2).*sqrt(1 + E.^2.*sin(2*theta).^2/(1 + E.*cos(theta).^2).^2 );
