function [han,hant,lon,lat] = plot_river_locations(file,g,color,rlist)
%
% [hansym,hanlab,lon,lat] = roms_plot_river_locations(file,g,[color],[rlist])
%
% Add symbols to an existing plot to show the locations of ROMS river 
% point sources.
%
% Symbols '>' and '<' indicate flow through a u-face in the positive,
%         negative direction, respectively, for the LuvSrc option
% Symbols '^' and 'v' indicate flow through a v-face in the positive,
%         negative direction, respectively, for the LuvSrc option
% Symbols 'o' indicate the file metadata states the sources are positioned
%         for the LwSrc option, and the symbols are plotted rho-point cell
%         centers. 
%
% Inputs:
%
% file     is a ROMS river forcing file
%          or a STRUCTURE containing the essential position information
%          generated using roms_get_river_source_locations.m
% g        is a ROMS grid structure computed using roms_get_grid.m
% color    (optional) is the symbol color
% rlist    (optional) vector of river numbers to plot (default is all)
%
% Outputs:
%
% hansym   is graphics handle to the plotted symbols (size number of rivers)
% hanlab   is graphics handle labels than number the rivers as in the file
% lon, lat are ROMS grid point coordinates of the river sources 
%

% John Wilkin; Sept 28 2014
%              Updated August 2016 to give more detail
%              Updated May 2020 to plot LwSrc sources.
%
% Copyright (c) 2021 - John L. Wilkin - jwilkin@rutgers.edu
% $Id: roms_plot_river_source_locations.m 595 2020-12-28 22:03:45Z wilkin $
%
% See also roms_get_river_source_locations.m

if 0
  % example usage
  % plot a checkerboard land/sea mask
  rmask = g.mask_rho;
  [mx,my] = meshgrid(1:size(rmask,2),1:size(rmask,1));
  cmask = rmask*2+mod(mod(mx,2)+mod(my,2),2)+1;
  pcolorjw(g.lon_rho,g.lat_rho,cmask); ff
  roms_plot_river_source_locations(file,g)
end

% Input FILE can be a ROMS river forcing file, or a structure containing
% the position information selected using roms_get_source_locations.m

if isstruct(file)
  xpos = file.xpos;
  ypos = file.epos;
  rdir = file.rdir;
  rsgn = file.rsgn;
else
  xpos = ncread(file,'river_Xposition');
  ypos = ncread(file,'river_Eposition');
  rdir = ncread(file,'river_direction');
  rsgn = sign(mean(ncread(file,'river_transport'),2,'omitnan'));
end

% check if this is a LwSrc formatted file
% I = ncinfo(file);
% if findstrinstruct(I.Variables,'Name','LwSrc')
%   LwSrc = ncread(file,'LwSrc');
%   if LwSrc
%     rdir = 2*ones(size(rdir));
%   end
% end

% get plot state

nextplotstatewas = get(gca,'nextplot');

% hold whatever is already plotted

set(gca,'nextplot','add')

% if nargin < 3
%   color = 'r';
% end

if nargin == 3
  if isnumeric(color)
    rlist = color;
    color = 'r';
  else
    rlist = 1:length(xpos);
  end
elseif nargin == 4
  if isnumeric(color)
    rlist = color;
  end
  if ischar(rlist)
    color = rlist;
  end
else
  color = 'r';
  rlist = 1:length(xpos);
end

% h = nan(size(xpos));
% ht = h;
% lon = h;
% lat = h;

kount = 0;
for r = rlist% 1:length(xpos)
  kount = kount+1;
  switch rdir(r)
    case 0
      lon(kount) = g.lon_u(xpos(r)+1, ypos(r));
      lat(kount) = g.lat_u(xpos(r)+1, ypos(r));
      switch rsgn(r)
        case 1
          sym = '>';
        case -1
          sym = '<';
        case 0
          sym = 's';
      end
      h(kount) = plot(lon(kount), lat(kount), [color sym]);
      % plot(lon(r),lat(r),[color '.']);
    case 1
      switch rsgn(r)
        case 1
          sym = '^';
        case -1
          sym = 'v';
        case 0
          sym = 's';
      end
      lon(kount) = g.lon_v(xpos(r), ypos(r)+1);
      lat(kount) = g.lat_v(xpos(r), ypos(r)+1);
      h(kount) = plot(lon(kount), lat(kount), [color sym]);
      % plot(lon(r),lat(r),[color 'x']);
    case 2
      sym = 'o';
      lon(kount) = g.lon_rho(xpos(r)+1, ypos(r)+1);
      lat(kount) = g.lat_rho(xpos(r)+1, ypos(r)+1);
      h(kount) = plot(lon(kount), lat(kount), [color sym]);
      % plot(lon(r),lat(r),[color 'x']);
      
    otherwise
      error('invalid river_direction')
  end
end

set(h,'MarkerSize',10)
set(h,'LineWidth',1)
% set(h,'MarkerFaceColor',color)

% label with numbers

plabel = true;
if plabel
  kount = 0;
  for r = rlist % 1:length(xpos)
    kount = kount+1;
    ht(kount) = text(lon(kount), lat(kount),['  ' int2str(r)],  ...
                     'fontsize', 12);
  end
end

% restore nextplotstate to what it was

set(gca,'nextplot',nextplotstatewas);

if nargout > 0
  han = h;
  if nargout > 1
    hant = ht;
  end
end
