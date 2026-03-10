function S=extract_section(G, field, Xpath, Ypath, varargin)

%
% EXTRACT_SECTION:  Extract a cross-section fro a 3D ROMS field
%
% S=extract_section(G, field, Xpath, Ypath, npath)
%
% It extracts a speficied cross-section (Xpath, Ypath) from a ROMS
% 3D field. The horizontally interpolated values are returned at
% all levels of the s-coordinate system. It is intended for plotting
% elsewhere.
%
% On Input:
%
%    Ginp        ROMS grid structure (struct array) with depth arrays so
%                3D variables can be processed (see get_roms_grid.m).
%
%    field       ROMS field to interpolate from (2D or 3D array).
%                  Avoid passing field time records. This function cannot
%                  be used for temporal interpolation.
%
%    Xpath       Longitude (degrees_east) or Cartesian (m) coordinates in
%                  the XI-direction (1D array).
%
%    Ypath       Latitude (degrees_north) or Cartesian (m) coordinates in
%                  the ETA-direction (1D array).
%
%    npath       Number of section points, OPTIONAL (default 100)
%
%                  If Xpath and Ypath has only two values, the extraction
%                  coordinates are computed as:
%
%                  x = linspace(Xpath(1), Xpath(2), npath)
%                  y = linspace(Ypath(1), Ypath(2), npath)
%
% On Output:
%
%    S           Extracted field cross-section (struct array)
%
%                  S.Xpath(:)     Extracted X-coordinates path
%                  S.Ypath(:)     Extracted Y-coordinates path
%                  S.dis(:,N)     Cross-section distance (km)
%                  S.depth(:,N)   Cross-section depth (m, negative)
%                  S.Xgrd(:,N)    Cross-section X-grid (degrees or km)
%                  S.Ygrd(:,N)    Cross-section Y-grid (degrees or km)
%                  S.h(:)         Cross-section bathymetry (m, negative)
%                  S.value(:,N)   Cross-section value

% git $Id$
%======================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                              %
%    Licensed under a MIT/X style license                              %
%    See License_ROMS.md                            Hernan G. Arango   %
%======================================================================%

% Initialize.

S = struct('Xpath'   , [], 'Ypath'   , [], 'depth'  , [],           ...
           'dis'     , [], 'Xgrd'    , [], 'Ygrd'   , [],           ...
           'h'       , [], 'value'   , []);

switch numel(varargin)
  case 0
    npath = 100;
  case 1
    npath = varargin{1};
end

if ~isstruct(G)
  error(' EXTRACT_SECTION: G must be a ROMS Grid structure');
end

if (ndims(field) == 3)
  [Im,Jm,Km] = size(field);
else
  disp(' EXTRACT_SECTION: Input FIELD must be a 3D array');
  error([' Number of dimensions = ', num2str(ndims(field))]);
end

if (length(Xpath) == 2)
  S.Xpath = linspace(Xpath(1), Xpath(2), npath);
else
  S.Xpath = Xpath;
end

if (length(Ypath) == 2)
  S.Ypath = linspace(Ypath(1), Ypath(2), npath);
else
  S.Ypath = Ypath;
end

% Interpolation method in scatteredInterpolant.

method = 'linear';

%-----------------------------------------------------------------------
% Determine C-grid type variable and coordinates.
%-----------------------------------------------------------------------

[Lp,Mp]=size(G.h);
L=Lp-1;
M=Mp-1;

if ((Im == Lp) && (Jm == Mp) && isfield(G,'mask_rho'))
  mask = G.mask_rho;
  if (G.spherical)
    Xgrd = G.lon_rho;
    Ygrd = G.lat_rho;
  else
    Xgrd = G.x_rho;
    Ygrd = G.y_rho;
  end
  Zgrd = G.z_r;
  h = G.h;

elseif ((Im == L) && (Jm == M) && isfield(G,'mask_psi'))
  mask = G.mask_psi;
  if (G.spherical)
    Xgrd = G.lon_psi;
    Ygrd = G.lat_psi;
  else
    Xgrd = G.x_psi;
    Ygrd = G.y_psi;
  end
  Zgrd = G.z_p;
  h  = 0.25 .* (G.h(1:L,1:M ) + G.h(2:Lp,1:M ) +                     ...
                G.h(1:L,2:Mp) + G.h(2:Lp,2:Mp));

elseif ((Im == L) && (Jm == Mp) && isfield(G,'mask_u'))
  mask=G.mask_u;
  if (G.spherical)
    Xgrd = G.lon_u;
    Ygrd = G.lat_u;
  else
    Xgrd = G.x_u;
    Ygrd = G.y_u;
  end
  Zgrd = G.z_u;
  h  = 0.5 .* (G.h(1:L,1:Mp) + G.h(2:Lp,1:Mp));

elseif ((Im == Lp) && (Jm == M) && isfield(G,'mask_v'))
  mask = G.mask_v;
  if (G.spherical)
    Xgrd = G.lon_v;
    Ygrd = G.lat_v;
  else
    Xgrd = G.x_v;
    Ygrd = G.y_v;
  end
  Zgrd = G.z_v;
  h  = 0.5 .* (G.h(1:Lp,1:M) + G.h(1:Lp,2:Mp));

end

%--------------------------------------------------------------------------
% Extract field cross-section along specified slice path.
%--------------------------------------------------------------------------

% Compute path distance.

if (G.spherical)
  dis = cumsum([0; sw_dist(S.Ypath(:), S.Xpath(:), 'km')]);
  S.dis = repmat(dis,[1 Km]);
end

% Set grid coordinates for cross-section.

S.Xgrd = repmat(S.Xpath, [1 Km]);
S.Ygrd = repmat(S.Ypath, [1 Km]);

if (rank(S.Xgrd) == 1)
  S.Xgrd = reshape(S.Xgrd, length(S.Xpath), Km);
end
if (rank(S.Ygrd) == 1)
  S.Ygrd = reshape(S.Ygrd, length(S.Ypath), Km);
end

% Use scatteredInterpolant - ROMS horizontal coordinates are not plaid.

Fi = scatteredInterpolant(Xgrd(:), Ygrd(:), h(:), method);

% Interpolate bathymery.

S.h = Fi(S.Xpath(:), S.Ypath(:));
S.h = -S.h;

% Interpolate depth and field level-by-level.

S.value = nan([length(S.Xpath) Km]);
S.depth = nan([length(S.Xpath) Km]);

for k=1:Km
  fk = squeeze(field(:,:,k));
  zk = squeeze(Zgrd(:,:,k));

  Fi.Values = fk(:);    S.value(:,k) = Fi(S.Xpath, S.Ypath);
  Fi.Values = zk(:);    S.depth(:,k) = Fi(S.Xpath, S.Ypath);
end

return
