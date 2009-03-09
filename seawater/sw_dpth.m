function DEPTHM = sw_dpth(P,LAT)

% SW_DPTH    Depth from pressure
%
% USAGE:  dpth = sw_dpth(P,LAT)
%
% DESCRIPTION:
%    Calculates depth in metres from pressure in dbars.
%
% INPUT:  (all must have same dimensions)
%   P   = Pressure    [db]
%   LAT = Latitude in decimal degress north [-90..+90]
%         (lat may have dimensions 1x1 or 1xn where P(mxn).
%
% OUTPUT:
%  dpth = depth [metres]
%
% AUTHOR:  Phil Morgan 92-04-06  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Unesco 1983. Algorithms for computation of fundamental properties of 
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.

% svn $Id$
%===========================================================================
% SW_DPTH   $Revision$  $Date$
%           Copyright (C) CSIRO, Phil Morgan 1992.
%=========================================================================

% CALLER:  general purpose
% CALLEE:  none

%-------------
% CHECK INPUTS
%-------------
[mP,nP] = size(P);
[mL,nL] = size(LAT);
if mL==1 & nL==1
  LAT = LAT*ones(size(P));
end %if  

if (mP~=mL) | (nP~=nL)              % P & LAT are not the same shape
     if (nP==nL) & (mL==1)          % LAT for each column of P
        LAT = LAT( ones(1,mP), : ); %     copy LATS down each column
                                    %     s.t. dim(P)==dim(LAT)
     else
        error('sw_depth.m:  Inputs arguments have wrong dimensions')
     end %if
end %if

Transpose = 0;
if mP == 1  % row vector
   P         =  P(:);
   LAT       =  LAT(:);
   Transpose = 1;
end %if

%-------------
% BEGIN
%-------------
% Eqn 25, p26.  Unesco 1983.

DEG2RAD = pi/180;
c1 = +9.72659;
c2 = -2.2512E-5;
c3 = +2.279E-10;
c4 = -1.82E-15;
gam_dash = 2.184e-6;

LAT = abs(LAT);
X   = sin(LAT*DEG2RAD);  % convert to radians
X   = X.*X;
bot_line = 9.780318*(1.0+(5.2788E-3+2.36E-5*X).*X) + gam_dash*0.5*P;
top_line = (((c4*P+c3).*P+c2).*P+c1).*P;
DEPTHM   = top_line./bot_line;

if Transpose
   DEPTHM = DEPTHM';
end %if

return
