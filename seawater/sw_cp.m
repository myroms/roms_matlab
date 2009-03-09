function cp = sw_cp(S,T,P)

% SW_CP      Heat Capacity (Cp) of sea water
%
% USAGE: cp = sw_cp(S,T,P)
%
% DESCRIPTION:
%    Heat Capacity of Sea Water using UNESCO 1983 polynomial.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (IPTS-68)]
%   P = pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%
% OUTPUT:
%   cp = Specific Heat Capacity  [J kg^-1 C^-1] 
% 
% AUTHOR:  Phil Morgan 93-04-20  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of 
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%

% svn $Id$
%=========================================================================
% SW_CP  $Revision$  $Date$
%         Copyright (C) CSIRO, Phil Morgan 1993.
%=========================================================================

% CALLER: general purpose
% CALLEE: none

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=3
   error('sw_cp.m: Must pass 3 parameters')
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

%------
% BEGIN
%------
P = P/10; % to convert db to Bar as used in Unesco routines

%------------
% eqn 26 p.32
%------------
c0 = 4217.4;
c1 =   -3.720283;
c2 =    0.1412855;
c3 =   -2.654387e-3;
c4 =    2.093236e-5;

a0 = -7.64357;
a1 =  0.1072763;
a2 = -1.38385e-3;

b0 =  0.1770383;
b1 = -4.07718e-3;
b2 =  5.148e-5;

Cpst0 =  c0 + c1.*T + c2.*T.^2 + c3.*T.^3 + c4.*T.^4 + ...
        (a0 + a1.*T + a2.*T.^2).*S + ...
	(b0 + b1.*T + b2.*T.^2).*S.*sqrt(S);
    
%------------
% eqn 28 p.33
%------------
a0 = -4.9592e-1;
a1 =  1.45747e-2;
a2 = -3.13885e-4;
a3 =  2.0357e-6;
a4 =  1.7168e-8;

b0 =  2.4931e-4;
b1 = -1.08645e-5;
b2 =  2.87533e-7;
b3 = -4.0027e-9;
b4 =  2.2956e-11;

c0 = -5.422e-8;
c1 =  2.6380e-9;
c2 = -6.5637e-11;
c3 =  6.136e-13;

del_Cp0t0 =  (a0 + a1.*T + a2.*T.^2 + a3.*T.^3 + a4.*T.^4).*P +    ...
	     (b0 + b1.*T + b2.*T.^2 + b3.*T.^3 + b4.*T.^4).*P.^2 + ...   
             (c0 + c1.*T + c2.*T.^2 + c3.*T.^3).*P.^3;

%------------	 
% eqn 29 p.34
%------------
d0 =  4.9247e-3;
d1 = -1.28315e-4;
d2 =  9.802e-7;
d3 =  2.5941e-8;
d4 = -2.9179e-10;

e0 = -1.2331e-4;
e1 = -1.517e-6;
e2 =  3.122e-8;

f0 = -2.9558e-6;
f1 =  1.17054e-7;
f2 = -2.3905e-9;
f3 =  1.8448e-11;

g0 =  9.971e-8;

h0 =  5.540e-10;
h1 = -1.7682e-11;
h2 =  3.513e-13;

j1 = -1.4300e-12;
S3_2  = S.*sqrt(S);

del_Cpstp = [(d0 + d1.*T + d2.*T.^2 + d3.*T.^3 + d4.*T.^4).*S + ...
             (e0 + e1.*T + e2.*T.^2).*S3_2].*P                + ...
	    [(f0 + f1.*T + f2.*T.^2 + f3.*T.^3).*S            + ...
	     g0.*S3_2].*P.^2                                  + ...
	     [(h0 + h1.*T + h2.*T.^2).*S                      + ...
	     j1.*T.*S3_2].*P.^3;
     

cp = Cpst0 + del_Cp0t0 + del_Cpstp;

if Transpose
   cp = cp';
end %if

return
