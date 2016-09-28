function t = Tsurf_Vlineardepth(c0, cp, X)
%tt = Tsurf_Vlineardepth(c0, cp, X)
%
%this function computes the arrival times from a surface source at
%positions X, based on a depth-dependent velocity model which is linear in
%depth. The source must be at X=0.
%
%input:
%   c0:     velocity at the surface (mm/us).
%   cp:     derivative of velocity with depth (/us).
%   X:      array of positions (mm).
%
%output:
%   t:      array of arrival times (us).


%compute ray parameter based on positions X:
pp = 2./sqrt(4*c0^2 + cp^2.*X.^2);

%compute arrival times:
t = (-2/cp)*log( c0*pp ./ (1 + sqrt(1-c0^2*pp.^2)) );