function [C,F]=grid_refine(Ginp,Gout,Gfactor,Imin,Imax,Jmin,Jmax)

%
% GRID_REFINE:  Creates a finer resolution Grid NetCDF file
%
% [C,F]=grid_refine(Ginp,Gout,Gfactor,Imin,Jmin,Imax,Jmax)
%
% Given a coarse resolution Grid NetCDF file (Ginp), this function
% creates a finer resolution grid in the region specified by the
% coarser grid coordinates (Imin,Jmin) and (Imax,Jmax). Notice that
% (Imin,Jmin), and (Imax,Jmax) indices are in terms of the PHI-points
% because it actually define the physical boundaries of the refined
% grid. We will add in the future the extra points needed for the
% contact zones.  The grid refinement coefficient is specified with
% Gfactor.
%
% On Input:
%
%    Ginp       Input  coaser Grid NetCDF file name (character string)
%    Gout       Output finner Grid NetCDF file name (character string)
%    Gfactor    Grid refinement factor (3,5,7,...)
%    Imin       Coarse grid lower-left  I-coordinate (PSI-point)
%    Imax       Coarse grid upper-right I-coordinate (PSI-point)
%    Imax       Coarse grid upper-right I-coordinate (PSI-point)
%    Jmin       Coarse grid lower-left  J-coordinate (PSI-point)
%    Jmax       Coarse grid upper-right J-coordinate (PSI-point)
%
% On Output:
%
%    C          Coarse resolution Grid structure
%    F          Fine   resolution Grid structure
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
  
% Set spherical switch.
  
spherical=nc_read(Ginp, 'spherical');
if (ischar(spherical)),
  if (spherical == 'T' || spherical == 't'),
    spherical=true;
  end
end

if (spherical),
  C.spherical = 1;
  F.spherical = 1;
else
  C.spherical = 0;
  F.spherical = 0;
end

% Get grid lengths.

C.xl = nc_read(Ginp, 'xl');
C.el = nc_read(Ginp, 'el');

F.xl = 0;
F.el = 0;

% Check grid coordinates and land/sea masking arrays available.

got.lon_rho = false;
got.lat_rho = false;
got.lon_psi = false;
got.lat_psi = false;
got.lon_u   = false;
got.lat_u   = false;
got.lon_v   = false;
got.lat_v   = false;

got.x_rho = false;
got.y_rho = false;
got.x_psi = false;
got.y_psi = false;
got.x_u   = false;
got.y_u   = false;
got.x_v   = false;
got.y_v   = false;

got.angle = false;

got.rmask = false;
got.pmask = false;
got.umask = false;
got.vmask = false;

[vname,nvars] = nc_vname(Ginp);

for n=1:nvars,
  name = deblank(vname(n,:));
  switch (name),
    case 'angle'
      got.angle = true;
    case 'lon_rho'
      got.lon_rho = true;
    case 'lat_rho'
      got.lat_rho = true;
    case 'lon_psi'
      got.lon_psi = true;
    case 'lat_psi'
      got.lat_psi = true;
    case 'lon_u'
      got.lon_u = true;
    case 'lat_u'
      got.lat_u = true;
    case 'lon_v'
      got.lon_v = true;
    case 'lat_v'
      got.lat_v = true;
    case 'x_rho'
      got.x_rho = true;
    case 'y_rho'
      got.y_rho = true;
    case 'x_psi'
      got.x_psi = true;
    case 'y_psi'
      got.y_psi = true;
    case 'x_u'
      got.x_u = true;
    case 'y_u'
      got.y_u = true;
    case 'x_v'
      got.x_v = true;
    case 'y_v'
      got.y_v = true;
    case 'mask_rho'
      got.rmask = true;
    case 'mask_psi'
      got.pmask = true;
    case 'mask_u'
      got.umask = true;
    case 'mask_v'
      got.vmask = true;
  end
end

% Set fields to process.

field_list = {'pm', 'pn', 'dmde', 'dndx', 'f', 'h'};

if (got.angle),
  field_list = [field_list, 'angle'];
end

if (got.x_rho && got.y_rho),
  field_list = [field_list, 'x_rho', 'y_rho'];
end

if (got.x_psi && got.y_psi),
  field_list = [field_list, 'x_psi', 'y_psi'];
end

if (got.x_u && got.y_u),
  field_list = [field_list, 'x_u', 'y_u'];
end

if (got.x_v && got.y_v),
  field_list = [field_list, 'x_v', 'y_v'];
end

if (spherical),
  field_list = [field_list, 'lon_rho', 'lat_rho',                     ...
                            'lon_psi', 'lat_psi',                     ...
                            'lon_u',   'lat_u',                       ...
                            'lon_v',   'lat_v'];
end

if (got.rmask || got.pmask || got.umask || got.vmask),
  field_list = [field_list, 'mask_rho', 'mask_psi',                   ...
                            'mask_u',   'mask_v'];
end

% Set coarser grid fractional coordinates.

[Lp,Mp] = size(nc_read(Ginp,'h'));      % RHO-points

L = Lp-1;  Lm = L-1;
M = Mp-1;  Mm = M-1;

[XrC,YrC] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);
[XpC,YpC] = meshgrid(1.0:1:M     , 1.0:1:L     );
[XuC,YuC] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );
[XvC,YvC] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);

C.Lr = Lp;
C.Mr = Mp;

% Check refinement region.

if (Imin >= Imax),
  error([' GRID_REFINE: Imin >= Imax,    ',                           ...
         ' Imin = ', num2str(Imin),                                   ...
         ' Imax = ', num2str(Imax)]);
end

if (Jmin >= Jmax),
  error([' GRID_REFINE: Jmin >= Jmax,    ',                           ...
         ' Jmin = ', num2str(Jmin),                                   ...
         ' Jmax = ', num2str(Jmax)]);
end

if (Imax > L),
  error([' GRID_REFINE: Imax > L,    ',                               ...
         ' Imax = ', num2str(Imax),                                   ...
         ' L = ', num2str(L)]);
end

if (Jmax > M),
  error([' GRID_REFINE: Jmax > M,    ',                               ...
         ' Jmax = ', num2str(Jmax),                                   ...
         ' M = ', num2str(M)]);
end

% Set finner grid fractional coordinates.

delta = 1.0/Gfactor;
fac0  = -0.5*delta;
fac1  = +0.5*delta;

IpF = (Imin:delta:Imax);                           % PSI-points
JpF = (Jmin:delta:Jmax);                           % PSI-points
IrF = (Imin+fac0:delta:Imax+fac1);                 % RHO-points
JrF = (Jmin+fac0:delta:Jmax+fac1);                 % RHO-points

[XrF,YrF] = meshgrid(JrF,IrF);                     % RHO-points
[XpF,YpF] = meshgrid(JpF,IpF);                     % PSI-points
[XuF,YuF] = meshgrid(JrF,IpF);                     % U-points
[XvF,YvF] = meshgrid(JpF,IrF);                     % V-points

F.Lr = 0;
F.Mr = 0;

%----------------------------------------------------------------------------
% Interpolate grid variables.
%----------------------------------------------------------------------------
%
% The default strategy is to interpolate the Cartesian coordinates
% (x_rho, y_rho, ...) in meters using the nondimensional fractional
% coordinates (XrC, YrC, ...).  Then, the other variables are
% interpolated in term of the Cartesian coordinates between
% coarse and finner grids.
%
% However, some application grid NetCDF files may not have the
% Cartesian coordinates (x_rho, t_rho, ...) variables. Then, the
% coarse to fine interpolation is done in terms of spherical
% coordinates (lon_rho, lat_rho, ...).
%
% The intepolation needs to be cubic guaranttee nonuniform grid
% spacing and curvilinear grids.

disp(' ');
disp(['Interpolating from coarse to fine ...']);

% Grid locations at RHO-points.

if (got.x_rho && got.y_rho),
  C.x_rho = nc_read(Ginp, 'x_rho');
  F.x_rho = interp2(XrC, YrC, C.x_rho,                                ...
                    XrF, YrF, 'cubic');

  C.y_rho = nc_read(Ginp, 'y_rho');
  F.y_rho = interp2(XrC, YrC, C.y_rho,                                ...
                    XrF, YrF, 'cubic');

  if (got.lon_rho && got.lat_rho),
    C.lon_rho = nc_read(Ginp, 'lon_rho');
    F.lon_rho = griddata(C.x_rho, C.y_rho, C.lon_rho,                 ...
                         F.x_rho, F.y_rho, 'cubic');

    C.lat_rho = nc_read(Ginp, 'lat_rho');
    F.lat_rho = griddata(C.x_rho, C.y_rho, C.lat_rho,                 ...
                         F.x_rho, F.y_rho, 'cubic');
  end
  
elseif (got.lon_rho && got.lat_rho),

  C.lon_rho = nc_read(Ginp, 'lon_rho');
  F.lon_rho = interp2(XrC, YrC, C.lon_rho,                            ...
                      XrF, YrF, 'cubic');

  C.lat_rho = nc_read(Ginp, 'lat_rho');
  F.lat_rho = interp2(XrC, YrC, C.lat_rho,                            ...
                      XrF, YrF, 'cubic');

  if (got.x_rho && got.x_rho),
    C.x_rho = nc_read(Ginp, 'x_rho');
    F.x_rho = griddata(C.lon_rho, C.lon_rho, C.x_rho,                 ...
                       F.lon_rho, F.lon_rho, 'cubic');

    C.y_rho = nc_read(Ginp, 'lat_rho');
    F.y_rho = griddata(C.lat_rho, C.lat_rho, C.y_rho,                 ...
                       F.lon_rho, F.lon_rho, 'cubic');
  end

else
  
  error([' GRID_REFINE: unable to find coordinates at RHO-points'])
  
end

% Grid locations at PSI-points.

if (got.x_psi && got.y_psi),
  C.x_psi = nc_read(Ginp, 'x_psi');
  F.x_psi = interp2(XpC, YpC, C.x_psi,                                ...
                    XpF, YpF, 'cubic');

  C.y_psi = nc_read(Ginp, 'y_psi');
  F.y_psi = interp2(XpC, YpC, C.y_psi,                                ...
                    XpF, YpF, 'cubic');

  if (got.lon_psi && got.lat_psi),
    C.lon_psi = nc_read(Ginp, 'lon_psi');
    F.lon_psi = griddata(C.x_psi, C.y_psi, C.lon_psi,                 ...
                         F.x_psi, F.y_psi, 'cubic');

    C.lat_psi = nc_read(Ginp, 'lat_psi');
    F.lat_psi = griddata(C.x_psi, C.y_psi, C.lat_psi,                 ...
                         F.x_psi, F.y_psi, 'cubic');
  end
  
elseif (got.lon_psi && got.lat_psi),

  C.lon_psi = nc_read(Ginp, 'lon_psi');
  F.lon_psi = interp2(XpC, YpC, C.lon_psi,                            ...
                      XpF, YpF, 'cubic');

  C.lat_psi = nc_read(Ginp, 'lat_psi');
  F.lat_psi = interp2(XpC, YpC, C.lat_psi,                            ...
                      XpF, YpF, 'cubic');

  if (got.x_psi && got.x_psi),
    C.x_psi = nc_read(Ginp, 'x_psi');
    F.x_psi = griddata(C.lon_psi, C.lon_psi, C.x_psi,                 ...
                       F.lon_psi, F.lon_psi, 'cubic');

    C.y_psi = nc_read(Ginp, 'lat_psi');
    F.y_psi = griddata(C.lat_psi, C.lat_psi, C.y_psi,                 ...
                       F.lon_psi, F.lon_psi, 'cubic');
  end

else
  
  error([' GRID_REFINE: unable to find coordinates at PSI-points'])

end

% Grid locations at U-points.

if (got.x_u && got.y_u),
  C.x_u = nc_read(Ginp, 'x_u');
  F.x_u = interp2(XuC, YuC, C.x_u,                                    ...
                  XuF, YuF, 'cubic');

  C.y_u = nc_read(Ginp, 'y_u');
  F.y_u = interp2(XuC, YuC, C.y_u,                                    ...
                  XuF, YuF, 'cubic');

  if (got.lon_u && got.lat_u),
    C.lon_u = nc_read(Ginp, 'lon_u');
    F.lon_u = griddata(C.x_u, C.y_u, C.lon_u,                         ...
                       F.x_u, F.y_u, 'cubic');

    C.lat_u = nc_read(Ginp, 'lat_u');
    F.lat_u = griddata(C.x_u, C.y_u, C.lat_u,                         ...
                       F.x_u, F.y_u, 'cubic');
  end
  
elseif (got.lon_u && got.lat_u),

  C.lon_u = nc_read(Ginp, 'lon_u');
  F.lon_u = interp2(XuC, YuC, C.lon_u,                                ...
                    XuF, YuF, 'cubic');

  C.lat_u = nc_read(Ginp, 'lat_u');
  F.lat_u = interp2(XuC, YuC, C.lat_u,                                ...
                    XuF, YuF, 'cubic');

  if (got.x_u && got.x_u),
    C.x_u = nc_read(Ginp, 'x_u');
    F.x_u = griddata(C.lon_u, C.lon_u, C.x_u,                         ...
                     F.lon_u, F.lon_u, 'cubic');

    C.y_u = nc_read(Ginp, 'lat_u');
    F.y_u = griddata(C.lat_u, C.lat_u, C.y_u,                         ...
                     F.lon_u, F.lon_u, 'cubic');
  end

else
  
  error([' GRID_REFINE: unable to find coordinates at U-points'])

end

% Grid locations at V-points.

if (got.x_v && got.y_v),
  C.x_v = nc_read(Ginp, 'x_v');
  F.x_v = interp2(XvC, YvC, C.x_v,                                    ...
                  XvF, YvF, 'cubic');

  C.y_v = nc_read(Ginp, 'y_v');
  F.y_v = interp2(XvC, YvC, C.y_v,                                    ...
                  XvF, YvF, 'cubic');

  if (got.lon_v && got.lat_v),
    C.lon_v = nc_read(Ginp, 'lon_v');
    F.lon_v = griddata(C.x_v, C.y_v, C.lon_v,                         ...
                       F.x_v, F.y_v, 'cubic');

    C.lat_v = nc_read(Ginp, 'lat_v');
    F.lat_v = griddata(C.x_v, C.y_v, C.lat_v,                         ...
                       F.x_v, F.y_v, 'cubic');
  end
  
elseif (got.lon_v && got.lat_v),

  C.lon_v = nc_read(Ginp, 'lon_v');
  F.lon_v = interp2(XvC, YvC, C.lon_v,                                ...
                    XvF, YvF, 'cubic');

  C.lat_v = nc_read(Ginp, 'lat_v');
  F.lat_v = interp2(XvC, YvC, C.lat_v,                                ...
                    XvF, YvF, 'cubic');

  if (got.x_v && got.x_v),
    C.x_v = nc_read(Ginp, 'x_v');
    F.x_v = griddata(C.lon_v, C.lon_v, C.x_v,                         ...
                     F.lon_v, F.lon_v, 'cubic');

    C.y_v = nc_read(Ginp, 'lat_v');
    F.y_v = griddata(C.lat_v, C.lat_v, C.y_v,                         ...
                     F.lon_v, F.lon_v, 'cubic');
  end

else
  
  error([' GRID_REFINE: unable to find coordinates at V-points'])
  
end

% Get grid lengths.

if (got.x_psi && got.y_psi),
  F.xl = max(F.x_psi(:)) - min(F.x_psi(:));
  F.el = max(F.y_psi(:)) - min(F.y_psi(:));
else
  F.xl = 0;
  F.el = 0;
end

% Other grid variables. The inverse metrics pm and pn cannot be
% interpolated. They need to be recomputed.

if (got.angle),
  C.angle = nc_read(Ginp, 'f');
  F.angle = griddata(C.x_rho, C.y_rho, C.angle,                       ...
                     F.x_rho, F.y_rho, 'cubic');
end

C.f     = nc_read(Ginp, 'f');
F.f     = griddata(C.x_rho, C.y_rho, C.f,                             ...
                   F.x_rho, F.y_rho, 'cubic');

C.pm    = nc_read(Ginp, 'pm');
C.pn    = nc_read(Ginp, 'pn');

C.dndx  = nc_read(Ginp, 'dndx');
C.dmde  = nc_read(Ginp, 'dmde');

F.dndx  = zeros(size(F.x_rho));
F.dmde  = zeros(size(F.x_rho));

% Land/sea masking: use nondimensional fractional coordinates.

if (got.rmask || got.pmask || got.umask || got.vmask),
  C.mask_rho = nc_read(Ginp, 'mask_rho');
  F.mask_rho = interp2(XrC, YrC, C.mask_rho,                          ...
                       XrF, YrF, 'nearest');

  [F.mask_u, F.mask_v, F.mask_psi]=uvp_masks(F.mask_rho);
end

% Bathymetry.

C.hraw = nc_read(Ginp, 'hraw', 1);
F.hraw = griddata(C.x_rho, C.y_rho, C.hraw,                           ...
                  F.x_rho, F.y_rho, 'cubic');

C.h    = nc_read(Ginp, 'h');
F.h    = griddata(C.x_rho, C.y_rho, C.h,                              ...
                  F.x_rho, F.y_rho, 'cubic');

% Recompute some variables.

deg2rad = pi / 180.0;
omega   = 2.0 * pi * 366.25 / (24.0 * 3600.0 * 365.25);

if (spherical),
  F.f_new = 2.0 * omega * sin(deg2rad * F.lat_rho);
end

disp(' ');
if (spherical),
  GreatCircle = true;
  disp(['Computing grid spacing: great circle distances']);
else
  GreatCircle = false;
  disp(['Computing grid spacing: Cartesian distances']);
end

[F.pm, F.pn, F.dndx, F.dmde]=grid_metrics(F, GreatCircle);

%----------------------------------------------------------------------------
% Create finner resolution Grid NetCDF file and write out data.
%----------------------------------------------------------------------------

% Set number of grid points.

[LpF,MpF] = size(F.x_rho);              % RHO-points

disp(' ');
disp(['Number of points:',                                            ...
      ' Coarse = ', num2str(Lp) ,' x ', num2str(Mp),                  ...
      ',  fine = ', num2str(LpF),' x ', num2str(MpF)]);

% Create ROMS Grid NetCDF file.

NewFile = true;

status = c_grid(LpF, MpF, Gout, NewFile);

% Set global attributes.

status = nc_attadd(Gout, 'parent_grid', Ginp);
status = nc_attadd(Gout, 'parent_Imin', int32(Imin));
status = nc_attadd(Gout, 'parent_Imax', int32(Imax));
status = nc_attadd(Gout, 'parent_Jmin', int32(Jmin));
status = nc_attadd(Gout, 'parent_Jmax', int32(Jmax));

status = nc_attadd(Gout, 'refine_factor', int32(Gfactor));

status = nc_attdel(Gout, 'history');
history = ['GRID file created using Matlab script: grid_refine, ', date_stamp];
status = nc_attadd(Gout, 'history', history);

% Write out fine resolution grid variables.

disp(['Writing finer grid variables into: ', Gout]);
disp(' ');

% Write out fine resolution grid variables.

status = nc_write (Gout, 'spherical', F.spherical);
status = nc_write (Gout, 'xl', F.xl);
status = nc_write (Gout, 'el', F.el);
status = nc_write (Gout, 'hraw', F.hraw, 1);

for value = field_list,
  field = char(value);
  status = nc_write (Gout, field, F.(field));
end,

return


