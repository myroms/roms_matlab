function u = fpt(p, doInverse, dx, dy)

% fpt -- Fast Poisson transform.
%  fpt(p) solves laplacian(u) = p, for matrix u,
%   assuming periodic boundary conditions.
%  fpt(p, doInverse) performs the inverse-transform
%   if "doInverse" is logically TRUE.  Default =
%   FALSE.
%  fpt(p, doInverse, dx, dy) uses sample intervals
%   dx and dy.  Defaults = 1 and 1, respectively.
%
% Also see: fps, fft2, ifft2.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1998 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================

% Reference: Press, et al., Numerical Recipes,
%    Cambridge University Press, 1986 and later.
 
% Version of 23-Oct-1998 09:02:58.
% Updated    10-Feb-2000 20:51:38.

if nargin < 1, help(mfilename), return, end
if nargin < 2, doInverse = 0; end
if nargin < 3, dx = 1; end
if nargin < 4, dy = 1; end

% Compute Fourier weights.

[m, n] = size(p);

i = (0:m-1).' * ones(1, n);
j = ones(m, 1) * (0:n-1);

% The simple formula for unit sample intervals is:

%  weights = 2 * (cos(2*pi*i/m) + cos(2*pi*j/n) - 2);

%  We modify it with the actual dx and dy values,
%  where dx runs horizontally and dy runs vertically.

dx2 = dx*dx;
dy2 = dy*dy;
weights = 2 * (cos(2*pi*i/m)/dy2 + cos(2*pi*j/n)/dx2 - 1/dy2 - 1/dx2);
weights(1, 1) = 1;

% Solve.

if ~any(doInverse)
	u = ifft2(fft2(p) ./ weights);
else
	u = ifft2(fft2(p) .* weights);
end

if isreal(p), u = real(u); end
