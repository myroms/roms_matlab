function tgriddata(nData, nGrid)

% tgriddata -- Compare "griddata" with "griddata1".
%  tgriddata(nData, nGrid) compares "griddata" with "griddata1",
%   for nData points (default = 100) and nGrid cells along each
%   edge of a square array (default = 10).
 
% svn $Id$
%=======================================================================
% Copyright (C) 2000 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 24-Jul-2000 14:01:04.
% Updated    24-Jul-2000 14:35:01.

if nargin < 1, help(mfilename), end

if nargin < 1, nData = 100; end
if nargin < 2, nGrid = 10; end

if ischar(nData), nData = eval(nData); end
if ischar(nGrid), nGrid = eval(nGrid); end

nData = round(nData);
nGrid = round(nGrid);

x = rand(nData, 1);
y = rand(size(x));
z = x + y;

g = linspace(0, 1, nGrid+1);

[xi, yi] = meshgrid(g, g);

tic

zi = griddata(x, y, z, xi, yi, 'linear');

f = find(isnan(zi));
if any(f)
	tri = delaunay(x, y);
	d = dsearch(x, y, tri, xi(f), yi(f));
	zi(f) = z(d);
end

t(1) = toc;

subplot(1, 2, 1)
hold off
surf(xi, yi, zi)
hold on
plot3(x, y, z, '*')
hold off
title('griddata')

tic

zi = griddata1(x, y, z, xi, yi, [inf 1 1]);

t(2) = toc;

subplot(1, 2, 2)
hold off
surf(xi, yi, zi)
hold on
plot3(x, y, z, '*')
hold off
title('griddata1')

try, zoomsafe on, catch, end

disp([' ## griddata  -- Elapsed time: ' num2str(t(1))])
disp([' ## griddata1 -- Elapsed time: ' num2str(t(2))])

ratio = int2str(round(t(2)/t(1)));
disp([' ## GRIDDATA is --> ' ratio ' <-- times faster than GRIDDATA1.'])

set(gcf, 'Name', [mfilename ' ' int2str(nData) ' ' int2str(nGrid)])
figure(gcf)
