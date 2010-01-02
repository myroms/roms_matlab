function [status]=c_climatology(S);

%
% C_CLIMATOLOGY:  Create a ROMS climatology NetCDF file
%
% [status]=c_climatology(S)
%
% This function creates a ROMS climatology NetCDF file using specified
% parameters in structure array, S.
%
% On Input:
%
%    S           Initial condidions creation parameters (structure array):
%
%                  S.ncname           NetCDF file name
%                  S.spherical        Spherical grid switch
%                  S.Vtransform       Vertical transformation equation
%                  S.Lm               Number of interior RHO-points in X
%                  S.Mm               Number of interior RHO-points in Y
%                  S.N                Number of vertical levels
%                  S.NT               Number of active and passive tracers
%
%                  S.def_zeta         Switch to define sea surfac height
%                  S.def_v2d          Switch to define 2D momentum
%                  S.def_v3d          Switch to define 3D momentum
%                  S.def_temp         Switch to define potential temperature
%                  S.def_salt         Switch to define salinity
%
%                  S.zeta_time        Number of sea surface height records
%                  S.v2d_time         Number of 2D momentum records
%                  S.v3d_time         Number of 3D momentum records
%                  S.temp_time        Number of potential temperature records
%                  S.salt_time        Number of salinity records
%
% On Output:
%
%    status      Error flag.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%----------------------------------------------------------------------------
%  Set some NetCDF parameters.
%----------------------------------------------------------------------------

[ncglobal ]=mexnc('parameter', 'nc_global');
[ncdouble ]=mexnc('parameter', 'nc_double');
[ncunlim  ]=mexnc('parameter', 'nc_unlimited');
[ncint    ]=mexnc('parameter', 'nc_int');
[ncfloat  ]=mexnc('parameter', 'nc_float');
[ncchar   ]=mexnc('parameter', 'nc_char');

%----------------------------------------------------------------------------
%  Get climatology creation parameters.
%----------------------------------------------------------------------------

if (isfield(S,'ncname')),
  ncname=S.ncname;
else,
  error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: ncname, ', ...
	  'in structure array S']);
end,

if (isfield(S,'spherical')),
  spherical=S.spherical;
else,
  spherical=0;
end,

if (isfield(S,'Vtransform')),
  Vtransform=S.Vtransform;
else,
  error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: Vtransform, ', ...
	  'in structure array S']);
end,

if (isfield(S,'Lm')),
  Lp=S.Lm+2;
else,
  error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: Lm, ', ...
          'in structure array S']);
end,

if (isfield(S,'Mm')),
  Mp=S.Mm+2;
else,
  error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: Mm, ', ...
	  'in structure array S']);
end,

if (isfield(S,'N')),
  N=S.N;
else,
  error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: N, ', ...
          'in structure array S']);
end,

if (isfield(S,'NT')),
  NT=S.NT;
else,
  error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: NT, ', ...
          'in structure array S']);
end,

if (isfield(S,'def_zeta')),
  if (S.def_zeta & ~isfield(S,'zeta_time')),
    error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: zeta_time, ', ...
            'in structure array S']);
  end,
else,   
  S.def_zeta=0;
  S.zeta_time=0;
end,

if (isfield(S,'def_v2d')),
  if (S.def_v2d & ~isfield(S,'v2d_time')),
    error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: v2d_time, ', ...
            'in structure array S']);
  end,
else,   
  S.def_v2d=0;
  S.v2d_time=0;
end,

if (isfield(S,'def_v3d')),
  if (S.def_v3d & ~isfield(S,'v3d_time')),
    error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: v3d_time, ', ...
            'in structure array S']);
  end,
else,   
  S.def_v3d=0;
  S.v3d_time=0;
end,

if (isfield(S,'def_temp')),
  if (S.def_temp & ~isfield(S,'temp_time')),
    error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: temp_time, ', ...
            'in structure array S']);
  end,
else,   
  S.def_temp=0;
  S.temp_time=0;
end,

if (isfield(S,'def_salt')),
  if (S.def_salt & ~isfield(S,'salt_time')),
    error([ 'C_CLIMATOLOGY - Cannot find dimension parameter: salt_time, ', ...
            'in structure array S']);
  end,
else,   
  S.def_salt=0;
  S.salt_time=0;
end,
  
%----------------------------------------------------------------------------
%  Set dimensions.
%----------------------------------------------------------------------------

Dname.xr = 'xi_rho';       Dsize.xr = Lp;
Dname.xu = 'xi_u';         Dsize.xu = Lp-1;
Dname.xv = 'xi_v';         Dsize.xv = Lp;
Dname.xp = 'xi_psi';       Dsize.xp = Lp-1;
Dname.yr = 'eta_rho';      Dsize.yr = Mp;
Dname.yu = 'eta_u';        Dsize.yu = Mp;
Dname.yv = 'eta_v';        Dsize.yv = Mp-1;
Dname.yp = 'eta_psi';      Dsize.yp = Mp-1;
Dname.Nr = 's_rho';        Dsize.Nr = N;
Dname.Nw = 's_w';          Dsize.Nw = N+1;
Dname.NT = 'tracer';       Dsize.NT = NT;

if (S.def_zeta),
 Dname.zeta_time = 'zeta_time';    Dsize.zeta_time = S.zeta_time;
end,

if (S.def_v2d),
 Dname.v2d_time  = 'v2d_time';     Dsize.v2d_time  = S.v2d_time;
end,

if (S.def_v3d),
 Dname.v2d_time  = 'v3d_time';     Dsize.v3d_time  = S.v3d_time;
end,

if (S.def_temp),
 Dname.temp_time = 'temp_time';    Dsize.temp_time = S.temp_time;
end,

if (S.def_salt),
 Dname.salt_time = 'salt_time';    Dsize.salt_time = S.temp_time;
end,

%----------------------------------------------------------------------------
%  Set Variables.
%----------------------------------------------------------------------------

%  Vertical grid variables.

Vname.Vtransform  = 'Vtransform';
Vname.Vstretching = 'Vstretching';
Vname.theta_s     = 'theta_s';
Vname.theta_b     = 'theta_b';
Vname.Tcline      = 'Tcline';
Vname.hc          = 'hc';
Vname.s_rho       = 's_rho';
Vname.s_w         = 's_w';
Vname.Cs_r        = 'Cs_r';
Vname.Cs_w        = 'Cs_w';

%  Horizontal grid variables.

Vname.spherical   = 'spherical';
Vname.h           = 'h';

if (spherical),
  Vname.rlon      = 'lon_rho';
  Vname.rlat      = 'lat_rho';
  Vname.ulon      = 'lon_u';
  Vname.ulat      = 'lat_u';
  Vname.vlon      = 'lon_v';
  Vname.vlat      = 'lat_v';
else,
  Vname.rx        = 'x_rho';
  Vname.ry        = 'y_rho';
  Vname.ux        = 'x_u';
  Vname.uy        = 'y_u';
  Vname.vx        = 'x_v';
  Vname.vy        = 'y_v';
end,

%  climatology variables.

Vname.zeta_time   = 'zeta_time';
Vname.v2d_time    = 'v2d_time';
Vname.v3d_time    = 'v3d_time';
Vname.temp_time   = 'temp_time';
Vname.salt_time   = 'salt_time';

Vname.zeta        = 'zeta';
Vname.ubar        = 'ubar';
Vname.vbar        = 'vbar';
Vname.u           = 'u';
Vname.v           = 'v';
Vname.temp        = 'temp';
Vname.salt        = 'salt';

%----------------------------------------------------------------------------
%  Create initial conditions NetCDF file.
%----------------------------------------------------------------------------

[ncid,status]=mexnc('create',ncname,'clobber');
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: CREATE - unable to create file: ', ncname]);
  return
end,

%----------------------------------------------------------------------------
%  Define dimensions.
%----------------------------------------------------------------------------

[did.xr,status]=mexnc('def_dim',ncid,Dname.xr,Dsize.xr); 
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.xr]);
  return
end,

[did.xu,status]=mexnc('def_dim',ncid,Dname.xu,Dsize.xu);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.xu]);
  return
end,

[did.xv,status]=mexnc('def_dim',ncid,Dname.xv,Dsize.xv);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.xv]);
  return
end,

[did.yr,status]=mexnc('def_dim',ncid,Dname.yr,Dsize.yr);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.yr]);
  return
end,

[did.yu,status]=mexnc('def_dim',ncid,Dname.yu,Dsize.yu);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.yu]);
  return
end,

[did.yv,status]=mexnc('def_dim',ncid,Dname.yv,Dsize.yv);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.yv]);
  return
end,

[did.Nr,status]=mexnc('def_dim',ncid,Dname.Nr,Dsize.Nr);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.Nr]);
  return
end,

[did.Nw,status]=mexnc('def_dim',ncid,Dname.Nw,Dsize.Nw);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.Nw]);
  return
end,

[did.NT,status]=mexnc('def_dim',ncid,Dname.NT,Dsize.NT);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ',Dname.NT]);
  return
end,

if (S.def_zeta),
  [did.zeta_time,status]=mexnc('def_dim',ncid,Dname.zeta_time,Dsize.zeta_time);
  if (status ~= 0),
    error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ', ...
            Dname.zeta_time]);
    return
  end,
end,

if (S.def_v2d),
  [did.v2d_time,status]=mexnc('def_dim',ncid,Dname.v2d_time,Dsize.v2d_time);
  if (status ~= 0),
    error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ', ...
	    Dname.v2d_time]);
    return
  end,
end,

if (S.def_v3d),
  [did.v3d_time,status]=mexnc('def_dim',ncid,Dname.v3d_time,Dsize.v3d_time);
  if (status ~= 0),
    error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ', ...
            Dname.v3d_time]);
    return
  end,
end,

if (S.def_temp),
  [did.temp_time,status]=mexnc('def_dim',ncid,Dname.temp_time,Dsize.temp_time);
  if (status ~= 0),
    error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ', ...
	    Dname.temp_time]);     
    return
  end,
end,

if (S.def_salt),
  [did.salt_time,status]=mexnc('def_dim',ncid,Dname.salt_time,Dsize.salt_time);
  if (status ~= 0),
    error([ 'C_CLIMATOLOGY: DEF_DIM - unable to define dimension: ', ...
	    Dname.salt_time]);     
    return
  end,
end,


%----------------------------------------------------------------------------
%  Create global attributes.
%----------------------------------------------------------------------------

type='CLIMATOLOGY file';
lstr=max(size(type));
[status]=mexnc('put_att_text',ncid,ncglobal,'type',ncchar,lstr,type);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: PUT_ATT_TEXT - unable to global attribure: type.']);
  return
end,

history=['Climatology file using Matlab script: c_climatology, ',date_stamp];
lstr=max(size(history));
[status]=mexnc('put_att_text',ncid,ncglobal,'history',ncchar,lstr,history);
if (status ~= 0),
  error([ 'C_CLIMATOLOGY: PUT_ATT_TEXT - unable to global attribure: history.']);
  return
end,

%----------------------------------------------------------------------------
%  Define configuration variables.
%----------------------------------------------------------------------------

% Define spherical switch.

Var.name      = Vname.spherical;
Var.type      = ncchar;
Var.dimid     = [];
Var.long      = 'grid type logical switch';
Var.flag_str  = ['T, ' ...
                 'F'];
Var.meaning   = ['spherical ' ...
                 'Cartesian'];
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

% Define vertical coordinate variables.

Var.name      = Vname.Vtransform;
Var.type      = ncint;
Var.dimid     = [];
Var.long      = 'vertical terrain-following transformation equation';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.Vstretching;
Var.type      = ncint;
Var.dimid     = [];
Var.long      = 'vertical terrain-following stretching function';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.theta_s;
Var.type      = ncdouble;
Var.dimid     = [];
Var.long      = 'S-coordinate surface control parameter';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.theta_b;
Var.type      = ncdouble;
Var.dimid     = [];
Var.long      = 'S-coordinate bottom control parameter';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.Tcline;
Var.type      = ncdouble;
Var.dimid     = [];
Var.long      = 'S-coordinate surface/bottom layer width';
Var.units     = 'meter';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.hc;
Var.type      = ncdouble;
Var.dimid     = [];
Var.long      = 'S-coordinate parameter, critical depth';
Var.units     = 'meter';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.s_rho;
Var.type      = ncdouble;
Var.dimid     = [did.Nr];
Var.long      = 'S-coordinate at RHO-points';
Var.min       = -1;
Var.max       = 0;
Var.positive  = 'up';
if (Vtransform == 1),
  Var.stdname = 'ocena_s_coordinate_g1';
elseif (Vtransform == 2),
  Var.stdname = 'ocena_s_coordinate_g2';
end,
Var.formula   = 's: s_rho C: Cs_r eta: zeta depth: h depth_c: hc';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.s_w;
Var.type      = ncdouble;
Var.dimid     = [did.Nw];
Var.long      = 'S-coordinate at W-points';
Var.min       = -1;
Var.max       = 0;
Var.positive  = 'up';
if (Vtransform == 1),
  Var.stdname = 'ocena_s_coordinate_g1';
elseif (Vtransform == 2),
  Var.stdname = 'ocena_s_coordinate_g2';
end,
Var.formula   = 's: s_w C: Cs_w eta: zeta depth: h depth_c: hc';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.Cs_r;
Var.type      = ncdouble;
Var.dimid     = [did.Nr];
Var.long      = 'S-coordinate stretching function at RHO-points';
Var.min       = -1;
Var.max       = 0;
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

Var.name      = Vname.Cs_w;
Var.type      = ncdouble;
Var.dimid     = [did.Nw];
Var.long      = 'S-coordinate stretching function at W-points';
Var.min       = -1;
Var.max       = 0;
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Define bathymetry.

Var.name      = Vname.h;
Var.type      = ncdouble;
Var.dimid     = [did.yr did.xr];
Var.long      = 'bathymetry at RHO-points';
Var.units     = 'meter';
if (spherical),
  Var.coord   = 'lon_rho lat_rho';
else,
  Var.coord   = 'x_rho y_rho';
end,
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Define horizontal grid variables.

if (spherical),
  Var.name      = Vname.rlon;
  Var.type      = ncdouble;
  Var.dimid     = [did.yr did.xr];
  Var.long      = 'longitude of RHO-points';
  Var.units     = 'degree_east';
  Var.stdname   = 'longitude';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.rlat;
  Var.type      = ncdouble;
  Var.dimid     = [did.yr did.xr];
  Var.long      = 'latitute of RHO-points';
  Var.units     = 'degree_north';
  Var.stdname   = 'latitude';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.ulon;
  Var.type      = ncdouble;
  Var.dimid     = [did.yu did.xu];
  Var.long      = 'longitude of U-points';
  Var.units     = 'degree_east';
  Var.stdname   = 'longitude';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.ulat;
  Var.type      = ncdouble;
  Var.dimid     = [did.yu did.xu];
  Var.long      = 'latitute of U-points';
  Var.units     = 'degree_north';
  Var.stdname   = 'latitude';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.vlon;
  Var.type      = ncdouble;
  Var.dimid     = [did.yv did.xv];
  Var.long      = 'longitude of V-points';
  Var.units     = 'degree_east';
  Var.stdname   = 'longitude';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.vlat;
  Var.type      = ncdouble;
  Var.dimid     = [did.yv did.xv];
  Var.long      = 'latitute of V-points';
  Var.units     = 'degree_north';
  Var.stdname   = 'latitude';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

else,

  Var.name      = Vname.rx;
  Var.type      = ncdouble;
  Var.dimid     = [did.yr did.xr];
  Var.long      = 'X-location of RHO-points';
  Var.units     = 'meter';
  [varid,status]=nc_vdef(ncid,Var);
  clear Var

  Var.name      = Vname.ry;
  Var.type      = ncdouble;
  Var.dimid     = [did.yr did.xr];
  Var.long      = 'Y-location of RHO-points';
  Var.units     = 'meter';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.ux;
  Var.type      = ncdouble;
  Var.dimid     = [did.yu did.xu];
  Var.long      = 'X-location of U-points';
  Var.units     = 'meter';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.uy;
  Var.type      = ncdouble;
  Var.dimid     = [did.yu did.xu];
  Var.long      = 'Y-location of U-points';
  Var.units     = 'meter';
  Var.field     = [Vname.uy,', scalar'];
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.vx;
  Var.type      = ncdouble;
  Var.dimid     = [did.yv did.xv];
  Var.long      = 'X-location of V-points';
  Var.units     = 'meter';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.vy;
  Var.type      = ncdouble;
  Var.dimid     = [did.yv did.xv];
  Var.long      = 'Y-location of V-points';
  Var.units     = 'meter';
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
  
end,

%----------------------------------------------------------------------------
%  Define climatology variables.
%----------------------------------------------------------------------------

%  Time variables.

if (S.def_zeta),
  Var.name      = Vname.zeta_time;
  Var.type      = ncdouble;
  Var.dimid     = [did.zeta_time];
  Var.long      = 'time for sea surface height climatology';
  Var.units     = 'day';
  Var.calendar  = '360.0 days in every year';
  Var.cycle     = 360.0;
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_v2d),
  Var.name      = Vname.v2d_time;
  Var.type      = ncdouble;
  Var.dimid     = [did.v2d_time];
  Var.long      = 'time for 2D momentum climatology';
  Var.units     = 'day';
  Var.calendar  = '360.0 days in every year';
  Var.cycle     = 360.0;
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_v3d),
  Var.name      = Vname.v3d_time;
  Var.type      = ncdouble;
  Var.dimid     = [did.v3d_time];
  Var.long      = 'time for 3D momentum climatology';
  Var.units     = 'day';
  Var.calendar  = '360.0 days in every year';
  Var.cycle     = 360.0;
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_temp),
  Var.name      = Vname.temp_time;
  Var.type      = ncdouble;
  Var.dimid     = [did.temp_time];
  Var.long      = 'time for potential temperature climatology';
  Var.units     = 'day';
  Var.calendar  = '360.0 days in every year';
  Var.cycle     = 360.0;
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_salt),
  Var.name      = Vname.salt_time;
  Var.type      = ncdouble;
  Var.dimid     = [did.salt_time];
  Var.long      = 'time for salinity climatology';
  Var.units     = 'day';
  Var.calendar  = '360.0 days in every year';
  Var.cycle     = 360.0;
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

%  Climatology fields.

if (S.def_zeta),
  Var.name      = Vname.zeta;
  Var.type      = ncdouble;
  Var.dimid     = [did.zeta_time did.yr did.xr];
  Var.long      = 'sea surface height climatology';
  Var.units     = 'meter';
  Var.time      = Vname.zeta_time;
  if (spherical),
    Var.coord   = strcat([Vname.rlon,' ',Vname.rlat,' ',Vname.zeta_time]); 
  else,
    Var.coord   = strcat([Vname.rx,' ',Vname.ry,' ',Vname.zeta_time]); 
  end,
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_v2d),
  Var.name      = Vname.ubar;
  Var.type      = ncdouble;
  Var.dimid     = [did.v2d_time did.yu did.xu];
  Var.long      = 'vertically integrated u-momentum component climatology';
  Var.units     = 'meter second-1';
  Var.time      = Vname.v2d_time;
  if (spherical),
    Var.coord   = strcat([Vname.ulon,' ',Vname.ulat,' ',Vname.v2d_time]);
  else,
    Var.coord   = strcat([Vname.ux,' ',Vname.uy,' ',Vname.v2d_time]);
  end,
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.vbar;
  Var.type      = ncdouble;
  Var.dimid     = [did.v2d_time did.yv did.xv];
  Var.long      = 'vertically integrated v-momentum component climatology';
  Var.units     = 'meter second-1';
  Var.time      = Vname.v2d_time;
  if (spherical),
    Var.coord   = strcat([Vname.vlon,' ',Vname.vlat,' ',Vname.v2d_time]); 
  else,
    Var.coord   = strcat([Vname.vx,' ',Vname.vy,' ',Vname.v2d_time]); 
  end,
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_v3d),
  Var.name      = Vname.u;
  Var.type      = ncdouble;
  Var.dimid     = [did.v3d_time did.Nr did.yu did.xu];
  Var.long      = 'u-momentum component climatology';
  Var.units     = 'meter second-1';
  Var.time      = Vname.v3d_time;
  if (spherical),
    Var.coord   = strcat([Vname.ulon,' ',Vname.ulat,' ',Vname.s_rho,' ', ...
                          Vname.v3d_time]); 
  else,
    Var.coord   = strcat([Vname.ux,' ',Vname.uy,' ',Vname.s_rho,' ', ...
	                  Vname.v3d_time]); 
  end,
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var

  Var.name      = Vname.v;
  Var.type      = ncdouble;
  Var.dimid     = [did.v3d_time did.Nr did.yv did.xv];
  Var.long      = 'v-momentum component climatology';
  Var.units     = 'meter second-1';
  Var.time      = Vname.v3d_time;
  if (spherical),
    Var.coord   = strcat([Vname.vlon,' ',Vname.vlat,' ',Vname.s_rho,' ', ...
	                  Vname.v3d_time]); 
  else,
    Var.coord   = strcat([Vname.vx,' ',Vname.vy,' ',Vname.s_rho,' ', ...
	                  Vname.v3d_time]); 
  end,
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_temp),
  Var.name      = Vname.temp;
  Var.type      = ncdouble;
  Var.dimid     = [did.temp_time did.Nr did.yr did.xr];
  Var.long      = 'potential temperature climatology';
  Var.units     = 'Celsius';
  Var.time      = Vname.temp_time;
  if (spherical),
    Var.coord   = strcat([Vname.rlon,' ',Vname.rlat,' ',Vname.s_rho,' ', ...
	                  Vname.temp_time]); 
  else,
    Var.coord   = strcat([Vname.rx,' ',Vname.ry,' ',Vname.s_rho,' ',
                          Vname.temp_time]); 
  end,
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

if (S.def_salt),
  Var.name      = Vname.salt;
  Var.type      = ncdouble;
  Var.dimid     = [did.salt_time did.Nr did.yr did.xr];
  Var.long      = 'salinity climatology';
  Var.time      = Vname.salt_time;
  if (spherical),
    Var.coord   = strcat([Vname.rlon,' ',Vname.rlat,' ',Vname.s_rho,' ', ...
	                  Vname.salt_time]); 
  else,
    Var.coord   = strcat([Vname.rx,' ',Vname.ry,' ',Vname.s_rho,' ', ...
	                  Vname.salt_time]); 
  end,
  [varid,status]=nc_vdef(ncid,Var);
  if (status ~= 0), return, end,
  clear Var
end,

%----------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%----------------------------------------------------------------------------

[status]=mexnc('enddef',ncid);
if (status == -1),
  error([ 'C_CLIMATOLOGY: ENDDEF - unable to leave definition mode.']);
  return
end,

[status]=mexnc('close',ncid);
if (status == -1),
  error([ 'C_CLIMATOLOGY: CLOSE - unable to close NetCDF file: ', ncname]);
  return
end,

return

