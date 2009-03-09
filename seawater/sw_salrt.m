function rt = sw_salrt(T)

% SW_SALRT   Conductivity ratio   rt(T)     = C(35,T,0)/C(35,15,0)
%
% USAGE:  rt = sw_salrt(T)
%
% DESCRIPTION:
%    Equation rt(T) = C(35,T,0)/C(35,15,0) used in calculating salinity.
%    UNESCO 1983 polynomial.
%
% INPUT: 
%   T = temperature [degree C (IPTS-68)]
%
% OUTPUT:
%   rt = conductivity ratio  [no units] 
% 
% AUTHOR:  Phil Morgan 93-04-17  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonoff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of 
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%

% svn $Id$
%=========================================================================
% SW_SALRT  $Revision$  $Date$
%           Copyright (C) CSIRO, Phil Morgan 1993.
%=========================================================================

% CALLER: sw_salt
% CALLEE: none

% rt = rt(T) = C(35,T,0)/C(35,15,0)
% Eqn (3) p.7 Unesco.

c0 =  0.6766097;
c1 =  2.00564e-2;
c2 =  1.104259e-4;
c3 = -6.9698e-7;
c4 =  1.0031e-9;

rt = c0 + (c1 + (c2 + (c3 + c4.*T).*T).*T).*T;

return
