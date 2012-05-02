function C = fine2coarse(Ginp,Gout,Gfactor,varargin)

%
% FINE2COARSE:  Creates a coarser resolution ROMS Grid NetCDF file
%
% C = fine2coarse(Ginp,Gout,Gfactor,Imin,Imax,Jmin,Jmax)
%
% Given a fine resolution Grid NetCDF file (Ginp), this function creates
% a coarser resolution grid in the region specified by the finer grid
% coordinates (Imin,Jmin) and (Imax,Jmax). Notice that (Imin,Jmin), and
% (Imax,Jmax) indices are in terms of the PSI-points because it actually
% defines the physical boundaries of the coarser grid. The grid coarseness
% coefficient is specified with Gfactor.
%
% If Imin, Imax, Jmin, and Jmax is not provided, this function extracts
% the largest coarse grid possible within the finer grid.  If these
% PSI-indices are provided make sure that you pick the correct value
% from the following set for consistency between coarse and fine grids:
%
%    Imin, Imax      Any value from    Gfactor : Gfactor :L-Gfactor
%    Jmin, Jmax      Any value from    Gfactor : Gfactor :M-Gfactor
%
% where L and M are the number of PSI-points in the I- and J-directions.
%
% On Input:
%
%    Ginp       Input  finer   Grid NetCDF file name (string)
%    Gout       Output coarser Grid NetCDF file name (string)
%    Gfactor    Grid coarseness factor (3,5,7,9,11,13,15,...)
%    Imin       Finer grid lower-left  I-coordinate (PSI-point)
%    Imax       Finer grid upper-right I-coordinate (PSI-point)
%    Jmin       Finer grid lower-left  J-coordinate (PSI-point)
%    Jmax       Finer grid upper-right J-coordinate (PSI-point)
%
% On Output:
%
%    C          Coaser resolution Grid structure
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%
  
% If applicable, get larger grid structure.

if (~isstruct(Ginp)),
  F = get_roms_grid(Ginp);
else
  F = Ginp;
end

% Set data sampling indices (PSI-points).

[Lp,Mp] = size(F.h);

L = Lp-1;
M = Mp-1;

% Set offset "half" factor which is used to extract boundary points
% outside of the PSI-points perimeter defined by Imin, Imax, Jmin,
% and Jmax.

half = floor(Gfactor-1)/2;

Imin = 1;                           % Default PSI-points range to apply
Imax = L-floor(mod(L,Gfactor)/2);   % coarseness to the full input grid
Jmin = 1;
Jmax = M-floor(mod(M,Gfactor)/2);

switch numel(varargin)
  case 1
    Imin = varargin{1};
  case 2
    Imin = varargin{1};
    Imax = varargin{2};
  case 3
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
  case 4
    Imin = varargin{1};
    Imax = varargin{2};
    Jmin = varargin{3};
    Jmax = varargin{4};
end

% Check coarseness factor.

legal = abs([3:2:27] - Gfactor) < 4*eps;
if (~any(legal)),
  error([' FINE2COARSE: illegal coarseness factor, Gfactor = ',         ...
         num2str(Gfactor)]);
end

% Check extraction region.

if (Imin >= Imax),
  error([' FINE2COARSE: Imin >= Imax,    ',                             ...
         ' Imin = ', num2str(Imin),                                     ...
         ' Imax = ', num2str(Imax)]);
end

if (Jmin >= Jmax),
  error([' FINE2COARSE: Jmin >= Jmax,    ',                             ...
         ' Jmin = ', num2str(Jmin),                                     ...
         ' Jmax = ', num2str(Jmax)]);
end

if (Imin < 1),
  error([' FINE2COARSE: Imin < 1,   ',                                  ...
         ' Imin = ', num2str(Imax)]);
end

if (Imax > L),
  error([' FINE2COARSE: Imax > L,    ',                                 ...
         ' Imax = ', num2str(Imax),                                     ...
         ' L = ', num2str(L)]);
end

if (Jmin < 1),
  error([' FINE2COARSE: Jmin < 1,   ',                                  ...
         ' Jmin = ', num2str(Imax)]);
end

if (Jmax > M),
  error([' FINE2COARSE: Jmax > M,    ',                                 ...
         ' Jmax = ', num2str(Jmax),                                     ...
         ' M = ', num2str(M)]);
end

% Set grid variables to process.

grd_vars = {'h', 'f', 'angle', 'pm', 'pn', 'dndx', 'dmde',              ...
            'x_rho', 'y_rho', 'x_psi', 'y_psi',                         ...
            'x_u', 'y_u', 'x_v', 'y_v'};

if (F.spherical),
  grd_vars = [grd_vars, 'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',     ...
                        'lon_u', 'lat_u', 'lon_v', 'lat_v'];
end

grd_vars = [grd_vars, 'mask_rho', 'mask_psi', 'mask_u', 'mask_v'];

grd_vars = [grd_vars, 'hraw'];               % needs to be last

%--------------------------------------------------------------------------
% Extract coaser grid.  Only valid for C-grid...
%--------------------------------------------------------------------------

% Set fine grid ROMS indices to extract. Notice the zero lower bound due
% to C-grid type conventions in ROMS.

IpF = 1:1:L;  JpF = 1:1:M;
IrF = 0:1:L;  JrF = 0:1:M;
IuF = 1:1:L;  JuF = 0:1:M;
IvF = 0:1:L;  JvF = 1:1:M;

% Set extraction ranges. The "half" offset value is to extract RHO-points
% boundaries (all four edges), U-points southern and northern boundary
% edges, and V-points western and eastern boundary edges. All these points
% are outside the fine grid perimeter defined by Imin, Imax, Jmin, and
% Jmax.

IminP = Imin;            ImaxP = Imax;
IminR = Imin-(half+1);   ImaxR = Imax+(half+1);
IminU = Imin;            ImaxU = Imax;
IminV = Imin-(half+1);   ImaxV = Imax+(half+1);

JminP = Jmin;            JmaxP = Jmax;
JminR = Jmin-(half+1);   JmaxR = Jmax+(half+1);
JminU = Jmin-(half+1);   JmaxU = Jmax+(half+1);
JminV = Jmin;            JmaxV = Jmax;

% Set indices to extract.  The "one" offset here is to account for the
% zero Fortran index that it is illegal in Matlab. These arrays elements
% are shifted by one in Matlab environment. The "fac" offset in the MOD
% function is to insure the correct shift due to the zero Fortran index
% in the MOD so we unity at the correct point, eq(mod(x,Gfactor)). This
% is correct!!!  You may test it if you like.  It is a good algebra
% exercise.

fac = half+1;         % use floor for integer arithmetic

IindexP = eq(mod(IpF-fac, Gfactor), 0) & ge(IpF,IminP)   & le(IpF,ImaxP);
JindexP = eq(mod(JpF-fac, Gfactor), 0) & ge(JpF,JminP)   & le(JpF,JmaxP);

IindexR = eq(mod(IrF    , Gfactor), 0) & ge(IrF,IminR-1) & le(IrF,ImaxR);
JindexR = eq(mod(JrF    , Gfactor), 0) & ge(JrF,JminR-1) & le(JrF,JmaxR);

IindexU = eq(mod(IuF-fac, Gfactor), 0) & ge(IuF,IminU  ) & le(IuF,ImaxU);
JindexU = eq(mod(JuF    , Gfactor), 0) & ge(JuF,JminU-1) & le(JuF,JmaxU);

IindexV = eq(mod(IvF    , Gfactor), 0) & ge(IvF,IminV-1) & le(IvF,ImaxV);
JindexV = eq(mod(JvF-fac, Gfactor), 0) & ge(JvF,JminV  ) & le(JvF,JmaxV);

% Initialize several subsomain structure paameters.

C.grid_name = Gout;
C.Lm = length(IpF(IindexP))-1;
C.Mm = length(JpF(JindexP))-1;
C.spherical = F.spherical;

% Extract grid variables.

for value = grd_vars,
  field = char(value);
  got.(field) = false;
  if (isfield(F,field)),
    switch (field)
      case {'lon_psi', 'lat_psi', 'mask_psi', 'x_psi', 'y_psi'}
        C.(field) = F.(field)(IindexP,JindexP);
        got.(field) = true;     
      case {'lon_u', 'lat_u', 'mask_u', 'x_u', 'y_u'}
        C.(field) = F.(field)(IindexU,JindexU);
        got.(field) = true;     
      case {'lon_v', 'lat_v', 'mask_v', 'x_v', 'y_v'}
        C.(field) = F.(field)(IindexV,JindexV);
          got.(field) = true;     
      case {'hraw'}     
        C.(field) = F.(field)(IindexR,JindexR,:);
        got.(field) = true;     
      otherwise
        C.(field) = F.(field)(IindexR,JindexR);
        got.(field) = true;     
    end
  end  
end

% Get grid lengths.

if (got.x_psi && got.y_psi),
  C.xl = max(C.x_psi(:)) - min(C.x_psi(:));
  C.el = max(C.y_psi(:)) - min(C.y_psi(:));
else
  C.xl = 0;
  C.el = 0;
end

%--------------------------------------------------------------------------
% Create subdomain Grid NetCDF file and write out data.
%--------------------------------------------------------------------------

% Set number of grid points.

[Lp,Mp] = size(C.h);                 % RHO-points

disp(' ');
disp(['Number of points:',                                              ...
      ' Fine   Grid = ', num2str(F.Lm+2) ,' x ', num2str(F.Mm+2),       ...
      ',  Coarse Grid = ', num2str(C.Lm+2),' x ', num2str(C.Mm+2)]);

% Create ROMS Grid NetCDF file.

NewFile = true;

status = c_grid(Lp, Mp, Gout, NewFile, C.spherical);
if (status ~= 0), return, end

% Set global attributes.

status = nc_attadd(Gout, 'parent_grid', Ginp);
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Imin', int32(Imin));
if (status ~= 0), return, end
  
status = nc_attadd(Gout, 'parent_Imax', int32(Imax));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Jmin', int32(Jmin));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'parent_Jmax', int32(Jmax));
if (status ~= 0), return, end

status = nc_attadd(Gout, 'coarse_factor', int32(Gfactor));
if (status ~= 0), return, end

status = nc_attdel(Gout, 'history');
if (status ~= 0), return, end

history = ['GRID file created using Matlab script: fine2coarse, ',      ...
           date_stamp];
status = nc_attadd(Gout, 'history', history);
if (status ~= 0), return, end

% Write out fine resolution grid variables.

disp(['Writing subdomain grid variables into: ', Gout]);
disp(' ');

status = nc_write (Gout, 'spherical', C.spherical);
if (status ~= 0), return, end

status = nc_write (Gout, 'xl', C.xl);
if (status ~= 0), return, end

status = nc_write (Gout, 'el', C.el);
if (status ~= 0), return, end

if (got.hraw),
  bath = size(C.hraw,3);
  for rec=1:bath,
    status = nc_write (Gout, 'hraw', C.hraw, rec);
    if (status ~= 0), return, end
  end
else
  status = nc_write (Gout, 'hraw', C.h, 1);
  if (status ~= 0), return, end
end

for value = 1:length(grd_vars)-1,
  field = char(grd_vars(value));
  if (got.(field)),
    status = nc_write (Gout, field, C.(field));
    if (status ~= 0), return, end
  end  
end,

%--------------------------------------------------------------------------
% Get full extracted grid structure.
%--------------------------------------------------------------------------

C = get_roms_grid(Gout);

return
