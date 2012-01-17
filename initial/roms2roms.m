function V = roms2roms(ncfile,P,T,Vname,Tindex,Rvector,varargin)

%
% ROMS2ROMS: Interpolates requested ROMS variable to specified ROMS grid
%
% V = roms2roms(ncfile,P,T,Vname,Tindex,Rvector,method,offset,RemoveNaN);
%
% This function interpolates requested 2D/3D variable between two ROMS
% application grids. The target grid must be inside of the parent grid.
%
% This function is intended for down-scaling or nesting applications.
% The horizontal/vertical coordinates for the parent and the target
% grids are specified with array structures 'P' and 'T', which are
% builded elsewhere using script 'get_roms_grid.m' for efficiency
% and functionality.  It uses 'TriScatteredInterp' for interpolating
% a 2D or 3D variable.
% 
% On Input:
%
%    ncfile        Parent NetCDF file/URL name (string) containing
%                    variable to process
%
%    P             Parent grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    T             Target grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    Vname         Field variable name to process (string)
%
%    Tindex        Time record index to process (scalar).
%
%    Rvector       Switch to Interpolate U- and V-points variables
%                    to RHO-points to facilitate rotation in
%                    in curvilinear grid applications somewhere
%                    else (logical)
%
%    method        Interpolation method in 'TriScatteredInterp'
%                    (string):
%
%                    'natural'     natural neighbor interpolation
%                    'linear'      linear interpolation (default)
%                    'nearest'     nearest-neighbor interpolation
%
%    offset        Number of extra points to used to sample the
%                    parent grid so is large enough to contain
%                    the target grid  (default 5)
%
%    RemoveNaN     Switch to remove NaN values from interpolated 
%                    variable with a second interpolation step
%                    using the nearest-neighbor method
%                    (default false)
%
% On Output:
%
%    V             Interpolated requested 2D/3D variable
%
 
% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set optional arguments.

method = 'linear';
offset = 5;
RemoveNaN = false;

switch numel(varargin)
  case 1
    method = varargin{1};
  case 2
    method = varargin{1};
    offset = varargin{2};
  case 3
    method = varargin{1};
    offset = varargin{2};
    RemoveNaN = varargin{3};
end

%  Initialize.

got.Mname = false;
got.Xname = false;
got.Yname = false;
got.Zname = false;

isr3d = false;
isw3d = false;
isvec = false;
recordless = true;

V = [];
Ncount = 0;

%  Get information about requested variable.
 
Info = nc_varinfo(ncfile,Vname);

nvdims = length(Info.Dimension);

%  Check variable dimensions and determine horizontal/vertical
%  coordinates and Land/Sea mask arrays.

if (nvdims > 0),
  for n=1:nvdims,
    dimnam = Info.Dimension{n};
    switch dimnam
      case 's_rho'
        isr3d = true;
      case 's_w'
        isw3d = true;
      case {'xi_rho','eta_rho'}
        Mname = 'mask_rho';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (P.spherical),
            Xname = 'lon_rho';
            Yname = 'lat_rho';
          else
            Xname = 'x_rho';
            Yname = 'y_rho';
          end
          Zname = 'z_r';
          got.Xname = true;
          got.Yname = true;
          got.Zname = true;
        end
      case {'xi_psi','eta_psi'}
        Mname = 'mask_psi';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (P.spherical),
            Xname = 'lon_psi';
            Yname = 'lat_psi';
          else
            Xname = 'x_psi';
            Yname = 'y_psi';
          end
          got.Xname = true;
          got.Yname = true;       
        end
      case {'xi_u','eta_u'}
        Mname = 'mask_u';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (P.spherical),
            Xname = 'lon_u';
            Yname = 'lat_u';
          else
            Xname = 'x_u';
            Yname = 'y_u';
          end 
          Zname = 'z_u';
          got.Xname = true;
          got.Yname = true;        
          got.Zname = true;
        end
        isvec = true;
     case {'xi_v','eta_v'}
        Mname = 'mask_v';
        got.Mname = true;
        if (~(got.Xname || got.Yname)),
          if (P.spherical),
            Xname = 'lon_v';
            Yname = 'lat_v';
          else
            Xname = 'x_v';
            Yname = 'y_v';
          end
          Zname = 'z_v';
          got.Xname = true;
          got.Yname = true;        
          got.Zname = true;
        end
        isvec = true;
      case 'ocean_time'
        recordless = false;    
    end
  end
  if (isw3d),
    Zname = 'z_w';
  end  
end

is3d = isr3d || isw3d;

%--------------------------------------------------------------------------
%  Get horizontal and vertical coordinates from parent and target grids.
%--------------------------------------------------------------------------

%  Parent grid.

if (isfield(P,Xname)),
  if (~isempty(P.(Xname)))
    XP = P.(Xname);
  else
    error([' ROMS2ROMS - field '', Xname, ''',                        ...
           ' is empty in parent grid structure: P']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Xname, ''',           ...
         ' in parent grid structure: P']);
end

if (isfield(P,Yname)),
  if (~isempty(P.(Yname)))
    YP = P.(Yname);
  else
    error([' ROMS2ROMS - field '', Yname, ''',                        ...
           ' is empty in parent grid structure: P']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Yname, ''',           ...
         ' in parent grid structure: P']);
end

if (is3d),
  if (isfield(P,Zname)),
    if (~isempty(P.(Zname)))
      ZP = P.(Zname);
    else
      error([' ROMS2ROMS - field '', Zname, ''',                      ...
             ' is empty in parent grid structure: P']);
    end
  else
    error([' ROMS2ROMS - unable to find field '', Zname, ''',         ...
           ' in parent grid structure: P']);
  end
end

if (isfield(P,Mname)),
  if (~isempty(P.(Mname)))
    Pmask = P.(Mname);
  else
    error([' ROMS2ROMS - field '', Mname, ''',                        ...
           ' is empty in parent grid structure: P']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Mname, ''',           ...
         ' in parent grid structure: P']);
end

%  Target grid.

if (isvec && Rvector),
  if (P.spherical),
    Xname = 'lon_rho';           % If requested, interpolate
    Yname = 'lat_rho';           % U- and V-points variables
  else                           % to RHO-points instead to
    Xname = 'x_rho';             % facilitate the curvilinear
    Yname = 'y_rho';             % rotation to true East and
  end                            % North.  This rotation is
  Mname = 'mask_rho';            % done somewhere else.
  if (is3d),
    Zname = 'z_r';
  end
end

if (isfield(T,Xname)),
  if (~isempty(T.(Xname)))  
    XT = T.(Xname);
  else
    error([' ROMS2ROMS - field '', Xname, ''',                        ...
           ' is empty in target grid structure: T']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Xname, ''',           ...
         ' in target grid structure: T']);
end

if (isfield(T,Yname)),
  if (~isempty(T.(Yname)))
    YT = T.(Yname);
  else
    error([' ROMS2ROMS - field '', Yname, ''',                        ...
           ' is empty in target grid structure: T']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Yname, ''',           ...
         ' in target grid structure: T']);
end

if (is3d),
  if (isfield(T,Zname)),
    if (~isempty(T.(Zname)))
      ZT = T.(Zname);
    else
      error([' ROMS2ROMS - field '', Zname, ''',                      ...
             ' is empty in target grid structure: T']);
    end
  else
    error([' ROMS2ROMS - unable to find field '', Zname, ''',         ...
           ' in target grid structure: T']);
  end
end
  
if (isfield(T,Mname)),
  if (~isempty(T.(Mname)))
    Tmask = T.(Mname);
  else
    error([' ROMS2ROMS - field '', Mname, ''',                        ...
           ' is empty in target grid structure: T']);
  end
else
  error([' ROMS2ROMS - unable to find field '', Mname, ''',           ...
         ' in target grid structure: T']);
end

%--------------------------------------------------------------------------
%  Read in requested variable from parent NetCDF file.
%--------------------------------------------------------------------------

ReplaceValue = NaN;
PreserveType = false;

VP = nc_read(ncfile,Vname,Tindex,ReplaceValue,PreserveType);

%--------------------------------------------------------------------------
%  Set parent grid sampling indices to accelerate the interpolation.
%  The horizontally sampled parent grid is large enough to contain
%  the target grid. The parameter 'offset' is used to add extra
%  points when computing the sampling indices (Istr:Iend,Jstr:Jend).
%  That is, the sampled grid is 'offset' points larger in all sides.
%  This is done to resolve well the interpolation near the boundaries
%  of the target grid.
%--------------------------------------------------------------------------

[Istr,Iend,Jstr,Jend] = sample_grid(XP,YP,XT,YT,offset);

%--------------------------------------------------------------------------
%  Interpolate requested variable to target grid.
%--------------------------------------------------------------------------

%  Determine if interpolating 2D or 3D variables.  Notice that it is
%  possible to interpolate recordless 2D or 3D variables.  That is,
%  variables with no time dimension.

if (recordless),
  ivdims = nvdims;
else
  ivdims = nvdims-1;
end

switch (ivdims),
 
 case 2

   [ImP,JmP]=size(XP);
   [ImT,JmT]=size(XT);

   disp(' ');
   disp(['Interpolating 2D variable: ', Vname,                        ...
          ' (', num2str(ImT), 'x', num2str(JmT),') from parent ',     ...
           '(', num2str(ImP), 'x', num2str(JmP),') ...']);
   disp(' ');

   x = XP(Istr:1:Iend,Jstr:1:Jend);
   y = YP(Istr:1:Iend,Jstr:1:Jend);
   v = VP(Istr:1:Iend,Jstr:1:Jend);

   x = x(:);
   y = y(:);
   v = v(:);

   Pind = find(Pmask(Istr:1:Iend,Jstr:1:Jend) < 0.5);    
   if (~isempty(Pind)),
     x(Pind) = [];                     % remove land points, if any
     y(Pind) = [];
     v(Pind) = [];
   end
   Pmin = min(v);
   Pmax = max(v);

   F = TriScatteredInterp(x,y,v,method);
   V = F(XT,YT);
   Tmin = min(V(:));
   Tmax = max(V(:));
    
   Tind = find(Tmask < 0.5);
   if (~isempty(Tind)),
     V(Tind) = 0;
   end,

%  If applicable, remove interpolated variable NaNs values with a
%  nearest neighbor interpolant.

   ind = find(isnan(V));

   if (~isempty(ind)),
     if (RemoveNaN),
       R = TriScatteredInterp(x,y,v,'nearest');
       V(ind) = R(XT(ind),YT(ind));
       Tmin = min(Tmin, min(V(ind)));
       Tmax = max(Tmax, max(V(ind)));

       ind = find(isnan(V));
       if (~isempty(ind)),
         Ncount = length(ind);
       end       
     else
       Ncount = length(ind);
     end
   end
   
   disp(['   Parent Min = ', sprintf('%12.5e',Pmin), '  ',            ...
           ' Parent Max = ', sprintf('%12.5e',Pmax)]);
   disp(['   Target Min = ', sprintf('%12.5e',Tmin), '  ',            ...
           ' Target Max = ', sprintf('%12.5e',Tmax), '  ',            ...
           ' Nan count = ',  num2str(Ncount)]);
   
 case 3

   [ImP,JmP,KmP]=size(ZP);
   [ImT,JmT,KmT]=size(ZT);

   disp(' ');
   disp(['Interpolating 3D variable: ', Vname,                        ...
          ' (', num2str(ImT), 'x', num2str(JmT), 'x',                 ...
	        num2str(KmT), ') from parent ',                       ...
           '(', num2str(ImP), 'x', num2str(JmP), 'x',                 ...
	        num2str(KmP), ') ...']);
   disp(' ');
  
   x = XP(Istr:1:Iend,Jstr:1:Jend);
   x = repmat(x,[1,1,KmP]); 
   y = YP(Istr:1:Iend,Jstr:1:Jend);
   y = repmat(y,[1,1,KmP]);
   z = ZP(Istr:1:Iend,Jstr:1:Jend,1:KmP);
   v = VP(Istr:1:Iend,Jstr:1:Jend,1:KmP);
  
   x = x(:);
   y = y(:);
   z = z(:);
   v = v(:);
  
   mask = Pmask(Istr:1:Iend,Jstr:1:Jend);
   Pind = find(repmat(mask,[1,1,KmP]) < 0.5);
   if (~isempty(Pind)),
     x(Pind) = [];                     % remove land points, if any
     y(Pind) = [];
     z(Pind) = [];
     v(Pind) = [];
   end
   Pmin = min(v);
   Pmax = max(v);

   F = TriScatteredInterp(x,y,z,v,method);
   X = repmat(XT,[1,1,KmT]);
   Y = repmat(YT,[1,1,KmT]);
   V = F(X,Y,ZT);
   Tmin = min(V(:));
   Tmax = max(V(:));
    
   Tind = find(repmat(Tmask,[1,1,KmT]) < 0.5);
   if (~isempty(Tind)),
     V(Tind) = 0;
   end,

%  If applicable, remove interpolated variable NaNs values with a
%  nearest neighbor interpolant.

   ind = find(isnan(V));

   if (~isempty(ind)),
     if (RemoveNaN),
       R = TriScatteredInterp(x,y,z,v,'nearest');
       V(ind) = R(X(ind),Y(ind), ZT(ind));
       Tmin = min(Tmin, min(V(ind)));
       Tmax = max(Tmax, max(V(ind)));

       ind = find(isnan(V));
       if (~isempty(ind)),
         Ncount = length(ind);
       end       
     else
       Ncount = length(ind);
     end
   end

   disp(['   Parent Min = ', sprintf('%12.5e',Pmin), '  ',            ...
           ' Parent Max = ', sprintf('%12.5e',Pmax)]);
   disp(['   Target Min = ', sprintf('%12.5e',Tmin), '  ',            ...
           ' Target Max = ', sprintf('%12.5e',Tmax), '  ',            ...
           ' Nan count = ',  num2str(Ncount)]);

end

return
