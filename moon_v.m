function [rmoon, rtasc, decl] = moon_v(jd)
% moon_v  Geocentric equatorial position of the Moon — vectorized.
%
%   jd    — Julian date, scalar or 1×N row vector
%   rmoon — N×3 position matrix, Earth radii (same unit as Vallado moon.m)
%   rtasc — 1×N right ascension, rad
%   decl  — 1×N declination, rad
%
% Algorithm: Vallado 2022, Alg 31.  All operations element-wise.

deg2rad = pi / 180.0;
twopi   = 2.0 * pi;

ttdb = (jd - 2451545.0) / 36525.0;   % 1×N

eclplong = 218.32  + 481267.8813 .*ttdb ...
    + 6.29 .*sin( (134.9 + 477198.85 .*ttdb).*deg2rad ) ...
    - 1.27 .*sin( (259.2 - 413335.38 .*ttdb).*deg2rad ) ...
    + 0.66 .*sin( (235.7 + 890534.23 .*ttdb).*deg2rad ) ...
    + 0.21 .*sin( (269.9 + 954397.70 .*ttdb).*deg2rad ) ...
    - 0.19 .*sin( (357.5 +  35999.05 .*ttdb).*deg2rad ) ...
    - 0.11 .*sin( (186.6 + 966404.05 .*ttdb).*deg2rad );    % deg

eclplat =   5.13 .*sin( ( 93.3 + 483202.03 .*ttdb).*deg2rad ) ...
    + 0.28 .*sin( (228.2 + 960400.87 .*ttdb).*deg2rad ) ...
    - 0.28 .*sin( (318.3 +   6003.18 .*ttdb).*deg2rad ) ...
    - 0.17 .*sin( (217.6 - 407332.20 .*ttdb).*deg2rad );    % deg

hzparal = 0.9508  + 0.0518 .*cos( (134.9 + 477198.85 .*ttdb).*deg2rad ) ...
    + 0.0095 .*cos( (259.2 - 413335.38 .*ttdb).*deg2rad ) ...
    + 0.0078 .*cos( (235.7 + 890534.23 .*ttdb).*deg2rad ) ...
    + 0.0028 .*cos( (269.9 + 954397.70 .*ttdb).*deg2rad );  % deg

eclplong = rem( eclplong.*deg2rad, twopi );   % rad
eclplat  = rem( eclplat .*deg2rad, twopi );
hzparal  = rem( hzparal .*deg2rad, twopi );

obliquity = (23.439291 - 0.0130042 .*ttdb) .* deg2rad;   % rad

l = cos(eclplat) .* cos(eclplong);
m = cos(obliquity).*cos(eclplat).*sin(eclplong) - sin(obliquity).*sin(eclplat);
n = sin(obliquity).*cos(eclplat).*sin(eclplong) + cos(obliquity).*sin(eclplat);

magr  = 1.0 ./ sin(hzparal);   % Earth radii

rmoon = [l; m; n]' .* magr';   % N×3
rtasc = atan2(m, l);            % 1×N
decl  = asin(n);                % 1×N

end
