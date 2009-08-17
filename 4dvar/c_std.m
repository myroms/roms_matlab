function [status]=c_std(S);

%
% C_STD:  Creates ROMS error covariance standard deviation NetCDF file
%
% [status]=c_std(S)
%
% This function creates ROMS 4DVar error covariance standard deviation
% NetCDF file using specified parameters in structure array, S.
%
% On Input:
%
%    S           Standard deviation creation parameters (structure array):
%
%                  S.ncname           NetCDF file name
%                  S.spherical        Spherical grid switch
%                  S.Vtransform       Vertical transformation equation
%                  S.Lm               Number of interior RHO-points in X
%                  S.Mm               Number of interior RHO-points in Y
%                  S.N                Number of vertical levels
%                  S.do_zeta          Switch to define free-surface
%                  S.do_ubar          Switch to define 2D U-velocity
%                  S.do_vbar          Switch to define 2D V-velocity
%                  S.do_u             Switch to define 3D U-velocity
%                  S.do_v             Switch to define 3D V-velocity
%                  S.do_temp          Switch to define temperature
%                  S.do_salt          Switch to define salinity
%
% On Output:
%
%    status      Error flag.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
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
%  Get error covariance standard deviation creation parameters.
%----------------------------------------------------------------------------

if (isfield(S,'ncname')),
  ncname=S.ncname;
else,
  error([ 'C_STD - Cannot find dimension parameter: ncname, ', ...
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
  error([ 'C_STD - Cannot find dimension parameter: Vtransform, ', ...
	  'in structure array S']);
end,

if (isfield(S,'Lm')),
  Lp=S.Lm+2;
else,
  error([ 'C_STD - Cannot find dimension parameter: Lm, ', ...
	  'in structure array S']);
end,

if (isfield(S,'Mm')),
  Mp=S.Mm+2;
else,
  error([ 'C_STD - Cannot find dimension parameter: Mm, ', ...
          'in structure array S']);
end,

if (isfield(S,'N')),
  N=S.N;
else,
  error([ 'C_STD - Cannot find dimension parameter: N, ', ...
          'in structure array S']);
end,

%----------------------------------------------------------------------------
%  Set dimensions.
%----------------------------------------------------------------------------

Dname.xr   = 'xi_rho';       Dsize.xr   = Lp;
Dname.xu   = 'xi_u';         Dsize.xu   = Lp-1;
Dname.xv   = 'xi_v';         Dsize.xv   = Lp;
Dname.xp   = 'xi_psi';       Dsize.xp   = Lp-1;
Dname.yr   = 'eta_rho';      Dsize.yr   = Mp;
Dname.yu   = 'eta_u';        Dsize.yu   = Mp;
Dname.yv   = 'eta_v';        Dsize.yv   = Mp-1;
Dname.yp   = 'eta_psi';      Dsize.yp   = Mp-1;
Dname.Nr   = 's_rho';        Dsize.Nr   = N;
Dname.Nw   = 's_w';          Dsize.Nw   = N+1;
Dname.time = 'ocean_time';   Dsize.time = ncunlim;

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

%  Initial conditions variables.

Vname.time        = 'ocean_time';
Vname.zeta        = 'zeta';
Vname.ubar        = 'ubar';
Vname.vbar        = 'vbar';
Vname.u           = 'u';
Vname.v           = 'v';
Vname.temp        = 'temp';
Vname.salt        = 'salt';

%----------------------------------------------------------------------------
%  Create standard deviation NetCDF file.
%----------------------------------------------------------------------------

[ncid,status]=mexnc('create',ncname,'clobber');
if (status ~= 0),
  error([ 'C_STD: CREATE - unable to create file: ', ncname]);
  return
end,

%----------------------------------------------------------------------------
%  Define dimensions.
%----------------------------------------------------------------------------

[did.xr,status]=mexnc('def_dim',ncid,Dname.xr,Dsize.xr); 
if (status ~= 0),
  error([ 'C_STD: ncdimdef - unable to define dimension: ',Dname.xr]);
  return
end,

[did.xu,status]=mexnc('def_dim',ncid,Dname.xu,Dsize.xu);
if (status ~= 0),
  error([ 'C_STD: DEF_DIM - unable to define dimension: ',Dname.xu]);
  return
end,

[did.xv,status]=mexnc('def_dim',ncid,Dname.xv,Dsize.xv);
if (status ~= 0),
  error([ 'C_STD: DEF_DIM - unable to define dimension: ',Dname.xv]);
  return
end,

[did.yr,status]=mexnc('def_dim',ncid,Dname.yr,Dsize.yr);
if (status ~= 0),
  error(['C_STD: DEF_DIM - unable to define dimension: ',Dname.yr]);
  return
end,

[did.yu,status]=mexnc('def_dim',ncid,Dname.yu,Dsize.yu);
if (status ~= 0),
  error([ 'C_STD: DEF_DIM - unable to define dimension: ',Dname.yu]);
  return
end,

[did.yv,status]=mexnc('def_dim',ncid,Dname.yv,Dsize.yv);
if (status ~= 0),
  error([ 'C_STD: DEF_DIM - unable to define dimension: ',Dname.yv]);
  return
end,

[did.Nr,status]=mexnc('def_dim',ncid,Dname.Nr,Dsize.Nr);
if (status ~= 0),
  error([ 'C_STD: DEF_DIM - unable to define dimension: ',Dname.Nr]);
  return
end,

[did.Nw,status]=mexnc('def_dim',ncid,Dname.Nw,Dsize.Nw);
if (status ~= 0),
  error([ 'C_STD: DEF_DIM - unable to define dimension: ',Dname.Nw]);
  return
end,

[did.time,status]=mexnc('def_dim',ncid,Dname.time,Dsize.time);
if (status ~= 0),
  error([ 'C_STD: DEF_DIM - unable to define dimension: ',Dname.time]);
  return
end,

%----------------------------------------------------------------------------
%  Create global attributes.
%----------------------------------------------------------------------------

type='ROMS/TOMS 4DVAR error covariance standard deviation';
lstr=max(size(type));
[status]=mexnc('PUT_ATT_TEXT',ncid,ncglobal,'type',ncchar,lstr,type);
if (status ~= 0),
  error([ 'C_STD: PUT_ATT_TEXT - unable to global attribure: type.']);
  return
end,

history=['Standard deviation file using Matlab script: c_std, ',date_stamp];
lstr=max(size(history));
[status]=mexnc('put_att_text',ncid,ncglobal,'history',ncchar,lstr,history);
if (status ~= 0),
  error([ 'C_STD: PUT_ATT_TEXT - unable to global attribure: history.']);
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
%  Define standard deviation variables.
%----------------------------------------------------------------------------

Var.name      = Vname.time;
Var.type      = ncdouble;
Var.dimid     = [did.time];
Var.long      = 'time since initialization';
Var.units     = 'seconds since 0001-01-01 00:00:00';
Var.calendar  = '360.0 days in every year';
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

if (isfield(S,'do_zeta')),
  if (S.do_zeta),
    Var.name      = Vname.zeta;
    Var.type      = ncdouble;
    Var.dimid     = [did.time did.yr did.xr];
    Var.long      = 'free-surface  standard deviation';
    Var.units     = 'meter';
    Var.time      = Vname.time;
    if (spherical),
      Var.coord   = strcat([Vname.rlon,' ',Vname.rlat,' ',Vname.time]); 
    else,
      Var.coord   = strcat([Vname.rx,' ',Vname.ry,' ',Vname.time]); 
    end,
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

if (isfield(S,'do_ubar')),
  if (S.do_ubar),
    Var.name      = Vname.ubar;
    Var.type      = ncdouble;
    Var.dimid     = [did.time did.yu did.xu];
    Var.long      = 'vertically integrated u-momentum component standard deviation';
    Var.units     = 'meter second-1';
    Var.time      = Vname.time;
    if (spherical),
      Var.coord   = strcat([Vname.ulon,' ',Vname.ulat,' ',Vname.time]); 
    else,
      Var.coord   = strcat([Vname.ux,' ',Vname.uy,' ',Vname.time]); 
    end,
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

if (isfield(S,'do_vbar')),
  if (S.do_vbar),
    Var.name      = Vname.vbar;
    Var.type      = ncdouble;
    Var.dimid     = [did.time did.yv did.xv];
    Var.long      = 'vertically integrated v-momentum component standard deviation';
    Var.units     = 'meter second-1';
    Var.time      = Vname.time;
    if (spherical),
      Var.coord   = strcat([Vname.vlon,' ',Vname.vlat,' ',Vname.time]); 
    else,
      Var.coord   = strcat([Vname.vx,' ',Vname.vy,' ',Vname.time]); 
    end,
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

if (isfield(S,'do_u')),
  if (S.do_u),
    Var.name      = Vname.u;
    Var.type      = ncdouble;
    Var.dimid     = [did.time did.Nr did.yu did.xu];
    Var.long      = 'u-momentum component standard deviation';
    Var.units     = 'meter second-1';
    Var.time      = Vname.time;
    if (spherical),
      Var.coord   = strcat([Vname.ulon,' ',Vname.ulat,' ',Vname.s_rho,' ',Vname.time]); 
    else,
      Var.coord   = strcat([Vname.ux,' ',Vname.uy,' ',Vname.s_rho,' ',Vname.time]); 
    end,
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,


if (isfield(S,'do_v')),
  if (S.do_v),
    Var.name      = Vname.v;
    Var.type      = ncdouble;
    Var.dimid     = [did.time did.Nr did.yv did.xv];
    Var.long      = 'v-momentum component standard deviation';
    Var.units     = 'meter second-1';
    Var.time      = Vname.time;
    if (spherical),
      Var.coord   = strcat([Vname.vlon,' ',Vname.vlat,' ',Vname.s_rho,' ',Vname.time]); 
    else,
      Var.coord   = strcat([Vname.vx,' ',Vname.vy,' ',Vname.s_rho,' ',Vname.time]); 
    end,
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

if (isfield(S,'do_temp')),
  if (S.do_temp),
    Var.name      = Vname.temp;
    Var.type      = ncdouble;
    Var.dimid     = [did.time did.Nr did.yr did.xr];
    Var.long      = 'potential temperature standard deviation';
    Var.units     = 'Celsius';
    Var.time      = Vname.time;
    if (spherical),
      Var.coord   = strcat([Vname.rlon,' ',Vname.rlat,' ',Vname.s_rho,' ',Vname.time]); 
    else,
      Var.coord   = strcat([Vname.rx,' ',Vname.ry,' ',Vname.s_rho,' ',Vname.time]); 
    end,
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

if (isfield(S,'do_salt')),
  if (S.do_salt),
    Var.name      = Vname.salt;
    Var.type      = ncdouble;
    Var.dimid     = [did.time did.Nr did.yr did.xr];
    Var.long      = 'salinity  standard deviation';
    Var.time      = Vname.time;
    if (spherical),
      Var.coord   = strcat([Vname.rlon,' ',Vname.rlat,' ',Vname.s_rho,' ',Vname.time]); 
    else,
      Var.coord   = strcat([Vname.rx,' ',Vname.ry,' ',Vname.s_rho,' ',Vname.time]); 
    end,
    [varid,status]=nc_vdef(ncid,Var);
    if (status ~= 0), return, end,
    clear Var
  end,
end,

%----------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%----------------------------------------------------------------------------

[status]=mexnc('enddef',ncid);
if (status == -1),
  error([ 'C_STD: ENDDEF - unable to leave definition mode.']);
  return
end,

[status]=mexnc('close',ncid);
if (status == -1),
  error([ 'C_STD: CLOSE - unable to close NetCDF file: ', ncname]);
  return
end,

return

