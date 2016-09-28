function phi = group_angle(theta,E)
%phi = GROUP_ANGLE(theta,E)
%
%This function computes the group angle as a function of the input phase
%angle and t anisotropy parameter, assuming the simple elliptical
%anisotropy model:
%
%  V_phase (theta)  = Vh * (1 + E*cos(theta)^2 )
%
% where theta is the phase angle.
%
%input:
%    theta:  phase angle
%    E:      anisotropy parameter
%
%output:
%    phi:    group angle

phi = atan(...
    (tan(theta) - E*sin(2*theta)./(1+E*cos(theta).^2))./(1+2*E*sin(theta).^2./(1+E*cos(theta).^2))...
    );