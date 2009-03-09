function u = fps(p, isSlope, dx, dy)

% fps -- Fast Poisson solver with boundary values.
%  fps(p) solves laplacian(u) = p for u, assuming the
%   boundary values are given along the perimeter of
%   p and the laplacians are given in the interior.
%  fps(p, isSlope) solves using the perimeter of p
%   as values of slope, if isSlope is logically TRUE.
%  fps(p, isSlope, dx, dy) uses sample intervals dx
%   and dy.  Defaults = 1 and 1, respectively.
%  fps([m n]) demonstrates itself with an n-by-n array
%   (default = [20 20]), with an off-center spike.
% Note: for the speediest Fourier transform, use a p
%  whose size is a composite of small factors, then
%  add 1, such as 33.
%
% Also see: fpt.
 
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
% Updated    10-Feb-2000 20:48:26.

if nargin < 1, help(mfilename), p = 'demo'; end
if nargin < 2, isSlope = 0; end
if nargin < 3, dx = 1; end
if nargin < 4, dy = 1; end

if ischar(isSlope), isSlope = eval(isSlope); end
isSlope == any(isSlope(:));

if isequal(p, 'demo'), p = 20; end
if ischar(p), p = eval(p); end
if length(p) == 1
	if nargin < 2
		p(2) = p;
	else
		p(2) = isSlope;
	end
end

if length(p) == 2
	theSize = max(p, 4);
	p = zeros(theSize);
	[m, n] = size(p);
	p(ceil(m/3), ceil(n/3)) = 1;
	feval(mfilename, p, isSlope, dx, dy);
	return
end

if nargout < 1, tic; end

% Fold the boundary into source terms.

if ~isSlope
	theFactor = -1;   % Boundary-value scheme.
else
	theFactor = +2;   % Boundary-slope scheme.
end

q = p;
for i = 1:2
	q = q.';
	q(2:end-1, 2:end-3:end-1) = ...
		q(2:end-1, 2:end-3:end-1) + ...
			theFactor * q(2:end-1, 1:end-1:end);
end

% Extract the interior.

q = q(2:end-1, 2:end-1);

% Symmetry: odd if 'value'; even if 'slope'.

if ~isSlope
	theSign = -1;   % Odd-symmetry, sine-transform scheme.
else
	theSign = +1;   % Even-symmetry, cosine-transform scheme.
end

[m, n] = size(q);
q = [zeros(m, 1) q zeros(m, 1) theSign*fliplr(q)];
q = [zeros(1, 2*n+2); q; zeros(1, 2*n+2); theSign*flipud(q)];

% Fast Poisson Transform.  The array q is now twice the size
%  of the original p, minus two elements in each direction.
%  Thus, to invoke a power-of-two Fourier transform, use
%  sizes that themselves are a power-of-two plus one.

res = fpt(q, 0, dx, dy);

% Retain relevant piece.

if ~isSlope
	result = p;
	result(2:end-1, 2:end-1) = res(2:m+1, 2:n+1);
else
	result = res(1:m+2, 1:n+2);
end

% Output or plot solution.

if nargout > 0
	u = result;
else
	disp([' ## Elapsed time: ' num2str(toc) ' s.'])
	subplot(2, 2, 1), surf(p), axis tight
	title('Laplacian Data'), xlabel('x'), ylabel('y'), zlabel('p')
	subplot(2, 2, 2), surf(result), axis tight
	title('Poisson Solution'), xlabel('x'), ylabel('y'), zlabel('u')
	subplot(2, 2, 3), plot(p), axis tight
	xlabel('x'), ylabel('p')
	subplot(2, 2, 4), plot(result), axis tight
	xlabel('x'), ylabel('u')
	figure(gcf)
	set(gcf, 'Name',[mfilename ' ' int2str(length(p))])
	error_norm = norm(4*del2(result(2:end-1, 2:end-1))-p(2:end-1, 2:end-1));
	disp([' ## Error norm: ' num2str(error_norm)])
end
