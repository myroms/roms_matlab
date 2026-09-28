function S=plot_rivers(Gname, Rname)

%
% PLOT_RIVERS:  Plots ROMS river point sources locations
%
% S = plot_rivers(Gname, Rname);
%
% This function plots the river point source locations on the
% discrete land/sea map.
%
% On Input:
%
%    Gname         ROMS Grid NetCDF file/URL name (string)
%              or, an existing ROMS grid structure (struct)
%
%    Rname         ROMS Rivers NetCDF filename (string)
%              or, and existing river location structure (struct)
%
% On Output:
%
%    S             River location structure (struct)
%
% Adapted from John Wilkin functions "roms_plot_mesh.m" and
% "roms_plot_river_source_locations.m".
%

% git $Id$
%======================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                              %
%    Licensed under a MIT/X style license           John L. Wilkin     %
%    See License_ROMS.md                            Hernan G. Arango   %
%======================================================================%

% Initialize.

S = struct('ncfile'           , [],                                  ...  
           'Nrivers'          , [],                                  ...
           'lon'              , [], 'lat'              , [],         ...
           'xpos'             , [], 'epos'             , [],         ...
           'rdir'             , [], 'rsgn'             , [],         ...
           'hancg'            , [], 'hancr'            , [],         ...
           'hansym'           , [], 'hanlab'           , []);

% Set ROMS Grid structure.

if (~isstruct(Gname))
  G = get_roms_grid(Gname);
else
  G = Gname;
end

% Set river locations.

if (~isstruct(Rname))
  S.ncfile = Rname;

  S.xpos = ncread(Rname, 'river_Xposition');
  S.epos = ncread(Rname, 'river_Eposition');
  S.rdir = ncread(Rname, 'river_direction');
else
  if (isfield(Rname, 'ncfile'))
    S.ncfile = Rname.ncfile;
  end
  S.xpos = Rname.xpos;
  S.epos = Rname.epos;
  S.rdir = Rname.rdir;
end

S.Nrivers = length(S.xpos);
S.rsgn = ones(size(S.xpos));

% Plot land/sea mask at RHO-points.

figure;

pcolorjw(G.lon_rho, G.lat_rho, G.mask_rho);
colormap(flipud(mpl_Set3));

% Plot discreate masked coastline. If available, plot geografical
% costal outine.

hold on;
S.hancg = plot_roms_mesh(G, 'coast', 'k');

a = axis;
if (isfield(G, 'lon_coast'))
  S.hancr = plot(G.lon_coast, G.lat_coast, 'r');
end
hold off;
axis(a);

% Plot river locations and enumerate index.

[S.hansym,S.hanlab,S.lon,S.lat] = plot_river_locations(S, G);
set(S.hansym,'color','b','markersize',5,'linewidth',2);
    
return