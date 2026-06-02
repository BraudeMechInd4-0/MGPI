function [rsun, rtasc, decl] = sun_v(jd)
% sun_v  Geocentric equatorial position of the Sun — vectorized.
%
%   jd   — Julian date, scalar or 1×N row vector
%   rsun — N×3 position matrix, AU (same unit as Vallado sun.m)
%   rtasc — 1×N right ascension, rad
%   decl  — 1×N declination, rad
%
% Algorithm: Vallado 2022, Alg 29.  Valid 1950–2050, ~0.01 deg accuracy.
% All operations element-wise; scalar if-branches replaced by logical indexing.

deg2rad = pi / 180.0;
twopi   = 2.0 * pi;

tut1 = (jd - 2451545.0) / 36525.0;   % 1×N

meanlong = rem( 280.460 + 36000.771285 .* tut1, 360.0 );   % deg

ttdb        = tut1;
meananomaly = rem( (357.528 + 35999.0509575 .* ttdb) .* deg2rad, twopi );  % rad
meananomaly(meananomaly < 0) = meananomaly(meananomaly < 0) + twopi;

eclplong = rem( meanlong + 1.914666471 .*sin(meananomaly) ...
                         + 0.019994643 .*sin(2.0.*meananomaly), 360.0 );    % deg

obliquity = (23.439291 - 0.0130042 .* ttdb) .* deg2rad;    % rad
eclplong  = eclplong .* deg2rad;                            % rad

magr = 1.000140612  - 0.016708617 .*cos(meananomaly) ...
                    - 0.000139589 .*cos(2.0.*meananomaly);  % AU

rsun = [magr.*cos(eclplong); ...
        magr.*cos(obliquity).*sin(eclplong); ...
        magr.*sin(obliquity).*sin(eclplong)]';              % N×3

rtasc = atan( cos(obliquity).*tan(eclplong) );              % 1×N

eclplong(eclplong < 0) = eclplong(eclplong < 0) + twopi;

idx = abs(eclplong - rtasc) > pi*0.5;
rtasc(idx) = rtasc(idx) + 0.5*pi .* round((eclplong(idx) - rtasc(idx)) ./ (0.5*pi));

decl = asin( sin(obliquity).*sin(eclplong) );               % 1×N

end
