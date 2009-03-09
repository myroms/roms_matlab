function ADTG = sw_adtg(S,T,P)

% SW_ADTG    Adiabatic temperature gradient
%
% adtg = sw_adtg(S,T,P)
%
% DESCRIPTION:
%    Calculates adiabatic temperature gradient as per UNESCO 1983 routines.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78) ]
%   T = temperature [degree C (IPTS-68)]
%   P = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   ADTG = adiabatic temperature gradient [degree_C/db]
%
% AUTHOR:  Phil Morgan 92-04-03  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonoff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of 
%    seawater. Unesco Tech. Pap. in Mar. Sci., No. 44, 53 pp.  Eqn.(31) p.39
%
%    Bryden, H. 1973.
%    "New Polynomials for thermal expansion, adiabatic temperature gradient
%    and potential temperature of sea water."
%    DEEP-SEA RES., 1973, Vol20,401-408.
%

% svn $Id$
%===========================================================================
% SW_ADTG   $Revision$  $Date$
%           Copyright (C) CSIRO, Phil Morgan  1992.
%===========================================================================

%-------------
% CHECK INPUTS
%-------------
if nargin ~= 3
   error('sw_adtg.m: Must pass 3 parameters ')
end %if

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
[mp,np] = size(P);

  
% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('check_stp: S & T must have same dimensions')
end %if

% CHECK OPTIONAL SHAPES FOR P
if     mp==1  & np==1      % P is a scalar.  Fill to size of S
   P = P(1)*ones(ms,ns);
elseif np==ns & mp==1      % P is row vector with same cols as S
   P = P( ones(1,ms), : ); %   Copy down each column.
elseif mp==ms & np==1      % P is column vector
   P = P( :, ones(1,ns) ); %   Copy across each row
elseif mp==ms & np==ns     % PR is a matrix size(S)
   % shape ok 
else
   error('check_stp: P has wrong dimensions')
end %if
[mp,np] = size(P);
 

  
% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0;
if mp == 1  % row vector
   P       =  P(:);
   T       =  T(:);
   S       =  S(:);   

   Transpose = 1;
end %if
%***check_stp

%-------------
% BEGIN
%-------------
a0 =  3.5803E-5;
a1 = +8.5258E-6;
a2 = -6.836E-8;
a3 =  6.6228E-10;

b0 = +1.8932E-6;
b1 = -4.2393E-8;

c0 = +1.8741E-8;
c1 = -6.7795E-10;
c2 = +8.733E-12;
c3 = -5.4481E-14;

d0 = -1.1351E-10;
d1 =  2.7759E-12;

e0 = -4.6206E-13;
e1 = +1.8676E-14;
e2 = -2.1687E-16;

ADTG =      a0 + (a1 + (a2 + a3.*T).*T).*T ... 
         + (b0 + b1.*T).*(S-35)  ...
	 + ( (c0 + (c1 + (c2 + c3.*T).*T).*T) + (d0 + d1.*T).*(S-35) ).*P ...
         + (  e0 + (e1 + e2.*T).*T ).*P.*P;

if Transpose
   ADTG = ADTG';
end %if

return
