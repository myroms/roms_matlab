function B = obc_roms2roms(ncfile,P,T,VarList,Tindex,boundary,varargin)

%
% OBC_ROMS2ROMS: Interpolates requested ROMS variable to specified ROMS grid
%
% B = obc_roms2roms(ncfile,P,T,VarList,Tindex,boundary, ...
%                   method,offset,RemoveNaN);
%
% This function interpolates lateral boundary conditions variables
% between two ROMS application grids. The target grid must be inside
% of the parent grid.
%
% This function is intended for down-scaling or nesting applications.
% The horizontal/vertical coordinates for the parent and the target
% grids are specified with array structures 'P' and 'T', which are
% builded elsewhere using script 'get_roms_grid.m' for efficiency
% and functionality. It uses 'TriScatteredInterp' for interpolating
% lateral boundary variables.
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
%    VarList       List of variables names to process (cell array)
%
%    Tindex        Time record index to process (scalar).
%
%    boundary      Lateral boundary condition switches of the grid
%                    edges to process (struct array)
%
%                    boundary.west        Western  edge
%                    boundary.east        Eastern  edge
%                    boundary.south       Southern edge
%                    boundary.north       Northern edge
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

%  Check if 'vector_rotation' field is available in the parent and
%  target Grid structures.Determine. Vector variables require special
%  processing.

if (~isfield(P,'vector_rotation'))
  if (isfield(P,'angle')),
    if (min(P.angle(:)) == 0 && max(P.angle(:)) == 0),
      P.vector_rotation = false;
    else
      P.vector_rotation = true;
    end
  else
    P.vector_rotation = false;
  end
end

if (~isfield(P,'vector_rotation'))
  if (isfield(P,'angle')),
    if (min(P.angle(:)) == 0 && max(P.angle(:)) == 0),
      P.vector_rotation = false;
    else
      P.vector_rotation = true;
    end
  else
    P.vector_rotation = false;
  end
end

Ucomponent = false;
Vcomponent = false;

%  Set report format for boundary variables.

lstr = 0;
for var = VarList,
  lstr = max(lstr, length(char(var)));
end
frmt = strcat('%',num2str(lstr+6),'s');

%==========================================================================
%  Process every field in the input variable list.
%==========================================================================

for var = VarList,

  Vname = char(var);

%  Initialize.

  got.Mname = false;
  got.Xname = false;
  got.Yname = false;
  got.Zname = false;

  isr3d = false;
  isw3d = false;
  isvec = false;

  Ncount = 0;

%  Get information about variable to process.

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
          Ucomponent = true;
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
          Vcomponent = true;
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
      error([' OBC_ROMS2ROMS - field '', Xname, ''',                  ...
             ' is empty in parent grid structure: P']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Xname, ''',     ...
           ' in parent grid structure: P']);
  end

  if (isfield(P,Yname)),
    if (~isempty(P.(Yname)))
      YP = P.(Yname);
    else
      error([' OBC_ROMS2ROMS - field '', Yname, ''',                  ...
             ' is empty in parent grid structure: P']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Yname, ''',     ...
           ' in parent grid structure: P']);
  end

  if (is3d),
    if (isfield(P,Zname)),
      if (~isempty(P.(Zname)))
        ZP = P.(Zname);
      else
        error([' OBC_ROMS2ROMS - field '', Zname, ''',                ...
               ' is empty in parent grid structure: P']);
      end
    else
      error([' OBC_ROMS2ROMS - unable to find field '', Zname, ''',   ...
             ' in parent grid structure: P']);
    end
  end

  if (isfield(P,Mname)),
    if (~isempty(P.(Mname)))
      Pmask = P.(Mname);
    else
      error([' OBC_ROMS2ROMS - field '', Mname, ''',                  ...
             ' is empty in parent grid structure: P']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Mname, ''',     ...
           ' in parent grid structure: P']);
  end

%  Target grid.

  if (isvec && P.vector_rotation),
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
      error([' OBC_ROMS2ROMS - field '', Xname, ''',                  ...
             ' is empty in target grid structure: T']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Xname, ''',     ...
           ' in target grid structure: T']);
  end

  if (isfield(T,Yname)),
    if (~isempty(T.(Yname)))
      YT = T.(Yname);
    else
      error([' OBC_ROMS2ROMS - field '', Yname, ''',                  ...
             ' is empty in target grid structure: T']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Yname, ''',     ...
           ' in target grid structure: T']);
  end

  if (is3d),
    if (isfield(T,Zname)),
      if (~isempty(T.(Zname)))
        ZT = T.(Zname);
      else
        error([' OBC_ROMS2ROMS - field '', Zname, ''',                ...
               ' is empty in target grid structure: T']);
      end
    else
      error([' OBC_ROMS2ROMS - unable to find field '', Zname, ''',   ...
             ' in target grid structure: T']);
    end
  end
  
  if (isfield(T,Mname)),
    if (~isempty(T.(Mname)))
      Tmask = T.(Mname);
    else
      error([' OBC_ROMS2ROMS - field '', Mname, ''',                  ...
             ' is empty in target grid structure: T']);
    end
  else
    error([' OBC_ROMS2ROMS - unable to find field '', Mname, ''',     ...
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
%  Set target grid lateral boundary conditions locations.  If parent or
%  target grids are curvilinear, extract two 
%--------------------------------------------------------------------------

  if (isvec && (P.vector_rotation || T.vector_rotation)),
    TX.west  = XT(1:2,:);
    TY.west  = YT(1:2,:);
    TM.west  = Tmask(1:2,:);

    TX.east  = XT(end-1:end,:);
    TY.east  = YT(end-1:end,:);
    TM.east  = Tmask(end-1:end,:);

    TX.south = XT(:,1:2);
    TY.south = YT(:,1:2);
    TM.south = Tmask(:,1:2);

    TX.north = XT(:,end-1:end);
    TY.north = YT(:,end-1:end);
    TM.north = Tmask(:,end-1:end);

    if (is3d),
      TZ.west  = ZT(1:2,:,:);
      TZ.east  = ZT(end-1:end,:,:);
      TZ.south = ZT(:,1:2,:);
      TZ.north = ZT(:,end-1:end,:);
    end
  else
    TX.west  = XT(1,:);
    TY.west  = YT(1,:);
    TM.west  = Tmask(1,:);

    TX.east  = XT(end,:);
    TY.east  = YT(end,:);
    TM.east  = Tmask(end,:);

    TX.south = XT(:,1);
    TY.south = YT(:,1);
    TM.south = Tmask(:,1);

    TX.north = XT(:,end);
    TY.north = YT(:,end);
    TM.north = Tmask(:,end);

    if (is3d),
      TZ.west  = ZT(1,:,:);
      TZ.east  = ZT(end,:,:);
      TZ.south = ZT(:,1,:);
      TZ.north = ZT(:,end,:);
    end
  end
  
%--------------------------------------------------------------------------
%  Interpolate lateral boundary conditions to target grid.
%--------------------------------------------------------------------------

%  Determine if processing 2D or 3D ROMS state variables.

  switch (nvdims-1),
 
    case 2

      [ImP,JmP]=size(XP);
      [ImT,JmT]=size(XT);

      disp(' ');
      disp(['Interpolating 2D OBC variable(s): ', Vname,              ...
            ' (', num2str(ImT), 'x', num2str(JmT),') from parent ',   ...
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
        x(Pind) = [];                  % remove land points, if any
        y(Pind) = [];
        v(Pind) = [];
      end
      Pmin = min(v);
      Pmax = max(v);

      F = TriScatteredInterp(x,y,v,method);

      for var = {'west','east','south','north'},

        edge = char(var);
        field = strcat(Vname,'_',edge);

        if (boundary.(edge)),
          B.(field) = F(TX.(edge), TY.(edge));
          Tmin = min(B.(field)(:));
          Tmax = max(B.(field)(:));
    
          Tind = find(TM.(edge) < 0.5);
          if (~isempty(Tind)),
            B.(field)(Tind) = 0;
          end

%  If applicable, remove interpolated variable NaNs values with a
%  nearest neighbor interpolant.

          ind = find(isnan(B.(field)));

          if (~isempty(ind)),
            if (RemoveNaN),
              R = TriScatteredInterp(x,y,v,'nearest');

              B.(field)(ind) = R(TX.(edge)(ind), TY.(field)(ind));
              Tmin = min(Tmin, min(B.(field)(ind)));
              Tmax = max(Tmax, max(B.(field)(ind)));

              ind = find(isnan(B.(field)));
              if (~isempty(ind)),
                Ncount = length(ind);
              end       
            else
              Ncount = length(ind);
            end
          end

          disp(['   ',sprintf(frmt,field),':  ',                      ...
                'Parent Min = ', sprintf('%12.5e',Pmin), '   ',       ...
                'Parent Max = ', sprintf('%12.5e',Pmax)]);
          disp(['   ',sprintf(frmt,field),':  ',                      ...
                'Target Min = ', sprintf('%12.5e',Tmin), '   ',       ...
                'Target Max = ', sprintf('%12.5e',Tmax), '   ',       ...
                'Nan count = ',  num2str(Ncount)]);
        end
      end
 
    case 3

      [ImP,JmP,KmP]=size(ZP);
      [ImT,JmT,KmT]=size(ZT);

      disp(' ');
      disp(['Interpolating 3D OBC variable(s): ', Vname,              ...
            ' (', num2str(ImT), 'x', num2str(JmT), 'x',               ...
                  num2str(KmT), ') from parent ',                     ...
             '(', num2str(ImP), 'x', num2str(JmP), 'x',               ...
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
        x(Pind) = [];                  % remove land points, if any
        y(Pind) = [];
        z(Pind) = [];
        v(Pind) = [];
      end
      Pmin = min(v);
      Pmax = max(v);

      F = TriScatteredInterp(x,y,z,v,method);

      for var = {'west','east','south','north'},

        edge = char(var);
        field = strcat(Vname,'_',edge);

        if (boundary.(edge)),
          X = repmat(TX.(edge),[1,1,KmT]);
          Y = repmat(TY.(edge),[1,1,KmT]);
   
          B.(field) = F(X, Y, TZ.(edge));
          Tmin = min(B.(field)(:));
          Tmax = max(B.(field)(:));
    
          Tind = find(repmat(TM.(edge),[1,1,KmT]) < 0.5);
          if (~isempty(Tind)),
            B.(field)(Tind) = 0;
          end

%  If applicable, remove interpolated variable NaNs values with a
%  nearest neighbor interpolant.

          ind = find(isnan(B.(field)));

          if (~isempty(ind)),
            if (RemoveNaN),
              R = TriScatteredInterp(x,y,z,v,'nearest');

              B.(field)(ind) = R(X(ind),Y(ind), TZ.(edge)(ind));
              Tmin = min(Tmin, min(B.(field)(ind)));
              Tmax = max(Tmax, max(B.(field)(ind)));

              ind = find(isnan(B.(field)));
              if (~isempty(ind)),
                Ncount = length(ind);
              end       
            else
              Ncount = length(ind);
            end
          end

          disp(['   ',sprintf(frmt,field),':  ',                      ...
                'Parent Min = ', sprintf('%12.5e',Pmin), '   ',       ...
                'Parent Max = ', sprintf('%12.5e',Pmax)]);
          disp(['   ',sprintf(frmt,field),':  ',                      ...
                'Target Min = ', sprintf('%12.5e',Tmin), '   ',       ...
                'Target Max = ', sprintf('%12.5e',Tmax), '   ',       ...
                'Nan count = ',  num2str(Ncount)]);
        end
      end
  end

%  Process next state variable in the list.

end

%--------------------------------------------------------------------------
%  Rotate vector components.
%--------------------------------------------------------------------------

%  Vector components that require rotation were interpolated at
%  RHO-points.  They are interpolated over two points adjacent to
%  the boundary edge to allow averaging to the appropriate C-grid
%  location.

if ((Ucomponent && Vcomponent) &&                                     ...
    (P.vector_rotation || T.vector_rotation)),
  B = rotate_vectors(B,P,T,VarList,boundary);
end

return

function B = rotate_vectors(B,P,T,VarList,boundary)

% This function rotates vector components to target application grid.
% In order to allow rotation, the vector components were interpolated
% at RHO-points. There is data available over the two points adjacent
% (columns/rows) to the target application grid boundary.  After the
% rotations are carried, the vector components are averaged to their
% respective C-grid locations.
%
% This function assumes that the 'parent_angle' field is in the
% Target Grid structure, T, and was interpolated from the Parent
% Grid somewhere else.
%
% On Input:
%
%    B             Lateral boundary conditions structure containing
%                    interpolated variables (struct array)
%
%    P             Parent grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    T             Target grid structure containing all horizontal
%                    and vertical variables (struct array)
%
%    VarList       List of variables names to process (cell array)
%
%    boundary      Lateral boundary condition switches of the grid
%                    edges to process (struct array)
%
%                    boundary.west        Western  edge
%                    boundary.east        Eastern  edge
%                    boundary.south       Southern edge
%                    boundary.north       Northern edge
%
% On Output:
%
%    Bnew          Update lateral boundary conditions structure
%                    (struct array)
%

%  Set target grid dimensions at RHO-points.

Lr = T.Lm + 2;
Mr = T.Mm + 2;
N  = T.N;

%--------------------------------------------------------------------------
%  Rotate 2D momentum vector components.
%--------------------------------------------------------------------------

if (max(strcmp(VarList,'ubar')) && max(strcmp(VarList,'vbar'))),

%  Parent Grid Orientation: rotate from (XI,ETA) coordinates to
%                           TRUE East/North.

  if (P.vector_rotation),
    angle.west  = T.parent_angle(1:2,:);
    angle.east  = T.parent_angle(end-1:end,:);
    angle.south = T.parent_angle(:,1:2);
    angle.north = T.parent_angle(:,end-1:end);

    for var = {'west','east','south','north'},
      edge = char(var);
      ufield = strcat('ubar','_',edge);
      vfield = strcat('vbar','_',edge);
      if (boundary.(edge)),
        if (isfield(B,ufield) && isfield(B,vfield)),
          ubar = B.(ufield);
          vbar = B.(vfield);

          B.(ufield) = ubar .* cos(angle.(edge)) -                    ...
                       vbar .* sin(angle.(edge));
          B.(vfield) = vbar .* cos(angle.(edge)) +                    ...
                       ubar .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,     ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Target Grid Orientation: rotate from TRUE East/North to
%                           (XI,ETA) coordinates.

  if (T.vector_rotation),
    angle.west  = T.angle(1:2,:);
    angle.east  = T.angle(end-1:end,:);
    angle.south = T.angle(:,1:2);
    angle.north = T.angle(:,end-1:end);

    for var = {'west','east','south','north'},
      edge = char(var);
      ufield = strcat('ubar','_',edge);
      vfield = strcat('vbar','_',edge);
      if (boundary.(edge)),
        if (isfield(B,ufield) && isfield(B,vfield)),
          ubar = B.(ufield);
          vbar = B.(vfield);

          B.(ufield) = ubar .* cos(angle.(edge)) +                    ...
                       vbar .* sin(angle.(edge));
          B.(vfield) = vbar .* cos(angle.(edge)) -                    ...
                       ubar .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,     ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Average vector components to U-point and V-points locations and
%  apply Lan/Sea mask.

  if (P.vector_rotation || T.vector_rotation),
    umask.west  = T.mask_u(1,:);
    umask.east  = T.mask_u(end,:);
    umask.south = T.mask_u(:,1);
    umask.north = T.mask_u(:,end);

    vmask.west  = T.mask_v(1,:);
    vmask.east  = T.mask_v(end,:);
    vmask.south = T.mask_v(:,1);
    vmask.north = T.mask_v(:,end);
    
    for var = {'west','east','south','north'},
      edge = char(var);
      ufield = strcat('ubar','_',edge);
      vfield = strcat('vbar','_',edge);
      if (boundary.(edge)),
        ubar = B.(ufield);
        vbar = B.(vfield);
        switch edge
          case {'west','east'}
            B.(ufield) = 0.5 .* (ubar(1,:) +                          ...
                                 ubar(2,:));
            B.(vfield) = 0.5 .* (vbar(1,1:Mr-1) +                     ...
                                 vbar(1,2:Mr  ));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
          case {'south','north'}
            B.(ufield) = 0.5.*(ubar(1:Lr-1,1)+                        ...
                               ubar(2:Lr  ,1));
            B.(vfield) = 0.5.*(vbar(:,1)+                             ...
                               vbar(:,2));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
        end
      end
    end
  end
end

%--------------------------------------------------------------------------
%  Rotate 3D momentum vector components.
%--------------------------------------------------------------------------

if (max(strcmp(VarList,'u')) && max(strcmp(VarList,'v'))),


%  Parent Grid Orientation: rotate from (XI,ETA) coordinates to
%                           TRUE East/North.

  if (P.vector_rotation),
    angle.west  = repmat(T.parent_angle(1:2,:), [1,1,N]);
    angle.east  = repmat(T.parent_angle(end-1:end,:), [1,1,N]);
    angle.south = repmat(T.parent_angle(:,1:2), [1,1,N]);
    angle.north = repmat(T.parent_angle(:,end-1:end), [1,1,N]);

    for var = {'west','east','south','north'},
      edge = char(var);
      ufield = strcat('u','_',edge);
      vfield = strcat('v','_',edge);
      if (boundary.(edge)),
        if (isfield(B,ufield) && isfield(B,vfield)),
          u = B.(ufield);
          v = B.(vfield);

          B.(ufield) = u .* cos(angle.(edge)) -                       ...
                       v .* sin(angle.(edge));
          B.(vfield) = v .* cos(angle.(edge)) +                       ...
                       u .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,     ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Target Grid Orientation: rotate from TRUE East/North to
%                           (XI,ETA) coordinates.

  if (T.vector_rotation),
    angle.west  = repmat(T.angle(1:2,:), [1,1,N]);
    angle.east  = repmat(T.angle(end-1:end,:), [1,1,N]);
    angle.south = repmat(T.angle(:,1:2), [1,1,N]);
    angle.north = repmat(T.angle(:,end-1:end), [1,1,N]);

    for var = {'west','east','south','north'},
      edge = char(var);
      ufield = strcat('u','_',edge);
      vfield = strcat('v','_',edge);
      if (boundary.(edge)),
        if (isfield(B,ufield) && isfield(B,vfield)),
          u = B.(ufield);
          v = B.(vfield);

          B.(ufield) = u .* cos(angle.(edge)) +                       ...
                       v .* sin(angle.(edge));
          B.(vfield) = v .* cos(angle.(edge)) -                       ...
                       u .* sin(angle.(edge));
        else
          error([' Unable to find fields: ',ufield,' or ',vfield,     ...
                 ' in structure array ''B''']);
        end
      end
    end
  end

%  Average vector components to U-point and V-points locations and
%  apply Land/Sea mask.

  if (P.vector_rotation || T.vector_rotation),
    umask.west  = squeeze(repmat(T.mask_u(1,:), [1,1,N]));
    umask.east  = squeeze(repmat(T.mask_u(end,:), [1,1,N]));
    umask.south = squeeze(repmat(T.mask_u(:,1), [1,1,N]));
    umask.north = squeeze(repmat(T.mask_u(:,end), [1,1,N]));

    vmask.west  = squeeze(repmat(T.mask_v(1,:), [1,1,N]));
    vmask.east  = squeeze(repmat(T.mask_v(end,:), [1,1,N]));
    vmask.south = squeeze(repmat(T.mask_v(:,1), [1,1,N]));
    vmask.north = squeeze(repmat(T.mask_v(:,end), [1,1,N]));
    
    for var = {'west','east','south','north'},
      edge = char(var);
      ufield = strcat('u','_',edge);
      vfield = strcat('v','_',edge);
      if (boundary.(edge)),
        u = B.(ufield);
        v = B.(vfield);
        switch edge
          case {'west','east'}
            B.(ufield) = squeeze(0.5 .* (u(1,:,:) +                   ...
                                         u(2,:,:)));
            B.(vfield) = squeeze(0.5 .* (v(1,1:Mr-1,:) +              ...
                                         v(1,2:Mr  ,:)));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
          case {'south','north'}
            B.(ufield) = squeeze(0.5.*(u(1:Lr-1,1,:)+                 ...
                                       u(2:Lr  ,1,:)));
            B.(vfield) = squeeze(0.5.*(v(:,1,:)+                      ...
                                       v(:,2,:)));

            B.(ufield) = B.(ufield) .* umask.(edge);
            B.(vfield) = B.(vfield) .* vmask.(edge);
        end
      end
    end
  end
end

return
