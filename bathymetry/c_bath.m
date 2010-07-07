function [status]=c_bath(Im, Jm, Bname);

%
% C_BATH:  Create a topography NetCDF file
%
% [status]=c_bath(Im,Jm,Bname)
%
% This function creates a topography NetCDF file.
%
% On Input:
%
%    Im          Number of RHO-points in the X-direction.
%    Jm          Number of RHO-points in the Y-direction.
%    Bname       Topography NetCDF file name.
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

%---------------------------------------------------------------------------
%  Get some NetCDF parameters.
%---------------------------------------------------------------------------

[ncglobal]=mexnc('parameter','nc_global');
[ncdouble]=mexnc('parameter','nc_double');
[ncunlim ]=mexnc('parameter','nc_unlimited');
[ncint   ]=mexnc('parameter','nc_int');
[ncfloat ]=mexnc('parameter','nc_float');
[ncchar  ]=mexnc('parameter','nc_char');

%---------------------------------------------------------------------------
%  Inquire dimensions from a existing NeTCDF file.
%---------------------------------------------------------------------------

Dname.lon ='lon';   Dsize.lon =Im;       Vname.lon ='lon';
Dname.lat ='lat';   Dsize.lat =Jm;       Vname.lat ='lat';
Dname.bath='bath';  Dsize.bath=ncunlim;  Vname.bath='hraw';

%---------------------------------------------------------------------------
%  Create topography NetCDF file.
%---------------------------------------------------------------------------

[ncid,status]=mexnc('nccreate',Bname,'nc_write');
if (ncid == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: nccreate - unable to create file: ', Bname]);
  return
end,

%---------------------------------------------------------------------------
%  Create global attribute(s).
%---------------------------------------------------------------------------

type='GRID file';
lstr=max(size(type));
[status]=mexnc('ncattput',ncid,ncglobal,'type',ncchar,lstr,type);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncattput - unable to global attribure: type.']);
  return
end,
history=['GRID file using Matlab script: c_bath, ', date_stamp];
lstr=max(size(history));
[status]=mexnc('ncattput',ncid,ncglobal,'history',ncchar,lstr,history);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncattput - unable to global attribure: history.']);
  return
end,

%---------------------------------------------------------------------------
%  Define dimensions.
%---------------------------------------------------------------------------

[did.lon]=mexnc('ncdimdef',ncid,Dname.lon,Dsize.lon);
if (did.lon == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncdimdef - unable to define dimension: ',Dname.lon]);
end,

[did.lat]=mexnc('ncdimdef',ncid,Dname.lat,Dsize.lat);
if (did.lat == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncdimdef - unable to define dimension: ',Dname.lat]);
end,

[did.bath]=mexnc('ncdimdef',ncid,Dname.bath,Dsize.bath);
if (did.bath == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncdimdef - unable to define dimension: ',Dname.bath]);
end,

%---------------------------------------------------------------------------
%  Define variables.
%---------------------------------------------------------------------------

% Define spherical switch.

Var.name          = 'spherical';
Var.type          = ncint;
Var.dimid         = [];
Var.long_name     = 'grid type logical switch';
Var.flag_values   = [0 1];
Var.flag_meanings = ['Cartesian', blanks(1), ...
                     'spherical'];
[varid,status]=nc_vdef(ncid,Var);
if (status ~= 0), return, end,
clear Var

%  Longitude.

Var.name          = Vname.lon;
Var.type          = ncdouble;
Var.dimid         = [did.lat did.lon];
Var.long_name     = 'longitude';
Var.units         = 'degree_east';
Var.standard_name = 'longitude';
[varid,status]=nc_vdef(ncid,Var);
clear Var

%  Latitude.

Var.name          = Vname.lat;
Var.type          = ncdouble;
Var.dimid         = [did.lat did.lon];
Var.long_name     = 'latitute';
Var.units         = 'degree_north';
Var.standard_name = 'latitude';
[varid,status]=nc_vdef(ncid,Var);
clear Var

%  Topography.

Var.name          = Vname.bath;
Var.type          = ncdouble;
Var.dimid         = [did.bath did.lat did.lon];
Var.long_name     = 'topography';
Var.units         = 'meter';
Var.coordinates   = strcat([Vname.lon,' ',Vname.lat]);
[varid,status]=nc_vdef(ncid,Var);
clear Var

%---------------------------------------------------------------------------
%  Leave definition mode and close NetCDF file.
%---------------------------------------------------------------------------

[status]=mexnc('ncendef',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncendef - unable to leave definition mode.']);
  return
end,

[status]=mexnc('ncclose',ncid);
if (status == -1),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['C_BATH: ncclose - unable to close GRID NetCDF file: ', Bname]);
  return
end,

return

