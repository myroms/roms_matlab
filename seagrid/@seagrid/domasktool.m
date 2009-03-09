function theResult = domasktool(self)

% seagrid/domasktool -- Draw the mask editing tool.
%  domasktool(self) draws the mask editing tool for self,
%   a "seagrid" object, superimposed on the current grid.
%   Boundary modification is suspended while the mask
%   editing tool is alive.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 15-Jun-1999 11:12:13.
% Updated    23-Apr-2001 16:00:02.

if nargout > 0, theResult = []; end
if nargin < 1, help(mfilename), return, end

theBathymetryFlag = psget(self, 'itsBathymetryFlag');
theMaskingFlag = psget(self, 'itsMaskingFlag');
theMaskTool = psget(self, 'itsMaskTool');

if ~any(theMaskTool)
	
	busy
	
	g = psget(self, 'itsSpacedGrids');
	x = g{1};
	y = g{2};
	
	theMask = psget(self, 'itsLand');

% Resize the mask, if needed.

	if ~isequal(size(theMask), size(x)-1)
		if isempty(theMask) | 1
			theMask = zeros(size(x)-1);
		else
			tic
			h = helpdlg('New grid size; Mask is now interpolated.', 'Please Note');
			drawnow
			theInterpFcn = psget(self, 'itsInterpFcn');
			theInterpMethod = psget(self, 'itsInterpMethod');
			before = size(theMask);
			after = size(x)-1;
			theMaskRange = [min(theMask(:)) max(theMask(:))];
			for i = 1:2
				count = before(i);
				bef = linspace(0, 1, count);
				bef = linspace(1/count, 1-1/count, count);
				count = after(i);
				aft = linspace(0, 1, count);
				aft = linspace(1/count, 1-1/count, count);
				theMask = feval(theInterpFcn, ...
						bef, theMask, ...
						aft, theInterpMethod).';
			end
			theMask = theMask > mean(theMaskRange);
			while ishandle(h) & toc < 3
				drawnow
			end
			delete(h)
		end
		theLand = ~~theMask;
		theWater = ~theLand;
		psset(self, 'itsMask', theWater)
		psset(self, 'itsWater', theWater)
		psset(self, 'itsLand', theLand)
	end

	theMask(end+1, end+1) = 0;
	
	z = [];
	c = [];
	
	bathy = psget(self, 'itsGriddedBathymetry');
	
	if theBathymetryFlag & ~isempty(bathy)
		bathy(end+1, :) = bathy(end, :);
		bathy(:, end+1) = bathy(:, end);
		z = bathy;
		if (0)
			theMask(isnan(z)) = 1;   % Mask the NaN depths.
		else
			z(isnan(z)) = 0;   % Set NaNs to zero depth.
		end
	end
	
	result = grid2mask(x, y, z, c, theMask);
	
	theAxes = get(result, 'Parent');
	set(theAxes, 'CLim', [-1 2])   % Why do this?
	f = findobj('Type', 'surface', 'Tag', 'grid2mask');
	set(f, 'ButtonDownFcn', 'psevent ButtonDownFcn', ...)
					'Tag', 'masktool', 'EraseMode', 'xor');
	
	h = get(result, 'UserData');
	set(h.strips, 'Tag', 'masktool')
	set(h.edit, 'Tag', 'masktool')
	
	theMaskTool = result;
	psset(self, 'itsMaskTool', theMaskTool)
	
	set(gcf, 'Pointer', 'crosshair')

else
	
	busy
	f = psget(self, 'itsMaskTool');
	if (0)   % Obsolete.
		h = get(f, 'UserData');
		theMask = zeros(size(h));
		[m, n] = size(theMask);
		for j = 1:n
			for i = 1:m
				switch get(h(i, j), 'Selected');
				case 'on'
					theMask(i, j) = 1;
				otherwise
				end
			end
		end
	end
	
	theMask = grid2mask('done');
	theMask = theMask(1:end-1, 1:end-1);
	
	theLand = ~~theMask;
	theWater = ~theLand;
	
	psset(self, 'itsMask', theMask)
	psset(self, 'itsWater', theWater)
	psset(self, 'itsLand', theLand)
	
	psset(self, 'itsMaskTool', [])
	refresh
	idle

end

if nargout > 0, theResult = self; end
