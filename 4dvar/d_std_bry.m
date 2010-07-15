%
%  D_STD_BRY:  Driver script to compute and write 4D-Var open boundary
%              conditions standard deviations.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  open boundary standard deviations which are used to convert modeled
%  error correlations to error covariances.
%
%  There are many ways to compute the open boundary conditions standard
%  deviations.  It depends on the source of the lateral boundary
%  conditions data.
% 
%  As an example, we are extracting the standard deviations from an
%  existing initial conditions standard deviation file. That is, a
%  full grid NetCDF file.  This not a good way to do it but it is
%  done here to show the structure of such computation.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Set output standard deviation NetCDF file.

STDfile = 'wc13_std_b.nc';

%  Set input grid file.

my_root = '/home/arango/ocean/toms/repository/test';

GRDfile = fullfile(my_root, 'WC13/Data', 'wc13_grd.nc');

%  Set input initial conditions standard deviation.

my_root = '/home/arango/ocean/toms/repository/test';

INIfile = fullfile(my_root, 'WC13/Data', 'wc13_std_i.nc');

%  Set boundary edges indices.

iwest  = 1;
isouth = 2;
ieast  = 3;
inorth = 4;

%---------------------------------------------------------------------------
%  Extract open boundary conditions standard deviations from initial
%  conditions standard deviation file.
%---------------------------------------------------------------------------

%  Set NetCDF file creation parameters.

[vnames,nvars] = nc_vname(INIfile);

S.curvilinear = false;
S.masking     = false;
S.Vtransform  = 1;
S.Vstretching = 1;

for n=1:nvars,
  name=deblank(vnames(n,:));
  switch name
    case 'spherical'
      S.spherical = nc_read(INIfile, 'spherical');
      if (ischar(S.spherical)),
        if (S.spherical == 'T' | S.spherical == 't');
          S.spherical = 1;
        else,
          S.spherical = 0;
        end,
      end,
    case 'Vtransform'
      S.Vtransform  = nc_read(INIfile, 'Vtransform');
    case 'Vstretching'
      S.Vstretching = nc_read(INIfile, 'Vstretching');
    case 'angle'
      S.angle = nc_read(INIfile, 'angle');
      S.curvilinear = true;
    case 'mask_rho'
      S.rmask = nc_read(INIfile, 'mask_rho');
      S.umask = nc_read(INIfile, 'mask_u');
      S.vmask = nc_read(INIfile, 'mask_v');
      S.masking = true;
  end,
end,

S.title = 'California Current System, 1/3 degree resolution (WC13)';

S.grd_file = GRDfile;

[Lr,Mr,Nr] = size(nc_read(INIfile, 'temp', 1));

S.Lm = Lr - 2;                      % number of interior RHO x-points
S.Mm = Mr - 2;                      % number of interior RHO y-points
S.N  = Nr;                          % number of vertical RHO levels

IorJ = max(Lr, Mr);                 % maximum number of lateral points

S.do_zeta = true;                   % free-surface
S.do_ubar = true;                   % vertically integrated u-momentum
S.do_vbar = true;                   % vertically integrated v-momentum
S.do_u    = true;                   % u-momentum
S.do_v    = true;                   % v-momentum
S.do_temp = true;                   % temperature
S.do_salt = true;                   % salinity

%  Read in grid.

S.rlon = nc_read(INIfile, 'lon_rho');
S.rlat = nc_read(INIfile, 'lat_rho');

S.ulon = nc_read(INIfile, 'lon_u');
S.ulat = nc_read(INIfile, 'lat_u');

S.vlon = nc_read(INIfile, 'lon_v');
S.vlat = nc_read(INIfile, 'lat_v');

%  Extract open boundary conditions standard deviations for initial
%  conditons.  This is just a template so user can see how the
%  boundary arrays are built. It is up to the user to compute the
%  appropriate values for a particalar application.

S.ncname = STDfile;

f = nc_read(INIfile, 'zeta', 1);                      % free-surface
[Im,Jm] = size(f);   S.zeta_obs = zeros(IorJ,4);
S.zeta_std(1:Jm,iwest ) = f(1  ,1:Jm);
S.zeta_std(1:Jm,ieast ) = f(end,1:Jm);
S.zeta_std(1:Im,isouth) = f(1:Im,1  );
S.zeta_std(1:Im,inorth) = f(1:Im,end);

f = nc_read(INIfile, 'ubar', 1);                      % 2D u-momentum
[Im,Jm] = size(f);   S.ubar_obs = zeros(IorJ,4);
S.ubar_std(1:Jm  ,iwest ) = f(1  ,1:Jm);
S.ubar_std(1:Jm  ,ieast ) = f(end,1:Jm);
S.ubar_std(2:Im+1,isouth) = f(1:Im,1  );
S.ubar_std(2:Im+1,inorth) = f(1:Im,end);

f = nc_read(INIfile, 'vbar', 1);                      % 2D v-momentum
[Im,Jm] = size(f);   S.vbar_std = zeros(IorJ,4);
S.vbar_std(2:Jm+1,iwest ) = f(1  ,1:Jm);
S.vbar_std(2:Jm+1,ieast ) = f(end,1:Jm);
S.vbar_std(1:Im  ,isouth) = f(1:Im,1  );
S.vbar_std(1:Im  ,inorth) = f(1:Im,end);

f = nc_read(INIfile, 'u', 1);                         % 3D u-momentum
[Im,Jm,Km] = size(f);   S.u_std = zeros(IorJ,Km,4);
S.u_std(1:Jm,1:Km  ,iwest ) = f(1  ,1:Jm,1:Km);
S.u_std(1:Jm,1:Km  ,ieast ) = f(end,1:Jm,1:Km);
S.u_std(2:Im+1,1:Km,isouth) = f(1:Im,1  ,1:Km);
S.u_std(2:Im+1,1:Km,inorth) = f(1:Im,end,1:Km);

f = nc_read(INIfile, 'v', 1);                         % 3D v-momentum
[Im,Jm,Km]=size(f);   S.v_std = zeros(IorJ,Km,4);
S.v_std(2:Jm+1,1:Km,iwest ) = f(1  ,1:Jm,1:Km);
S.v_std(2:Jm+1,1:Km,ieast ) = f(end,1:Jm,1:Km);
S.v_std(1:Im,1:Km  ,isouth) = f(1:Im,1  ,1:Km);
S.v_std(1:Im,1:Km  ,inorth) = f(1:Im,end,1:Km);

f = nc_read(INIfile, 'temp', 1);                      % temperature
[Im,Jm,Km]=size(f);   S.temp_std = zeros(IorJ,Km,4);
S.temp_std(1:Jm,1:Km,iwest ) = f(1  ,1:Jm,1:Km);
S.temp_std(1:Jm,1:Km,ieast ) = f(end,1:Jm,1:Km);
S.temp_std(1:Im,1:Km,isouth) = f(1:Im,1  ,1:Km);
S.temp_std(1:Im,1:Km,inorth) = f(1:Im,end,1:Km);

f = nc_read(INIfile, 'salt', 1);                      % salinity
[Im,Jm,Km] = size(f);   S.salt_std = zeros(IorJ,Km,4);
S.salt_std(1:Jm,1:Km,iwest ) = f(1  ,1:Jm,1:Km);
S.salt_std(1:Jm,1:Km,ieast ) = f(end,1:Jm,1:Km);
S.salt_std(1:Im,1:Km,isouth) = f(1:Im,1  ,1:Km);
S.salt_std(1:Im,1:Km,inorth) = f(1:Im,end,1:Km);

%---------------------------------------------------------------------------
%  Write out standard deviation fields.
%---------------------------------------------------------------------------

disp(' ');
disp([ 'Writting out standard deviation, file = ', S.ncname]);
disp(' ');

%  Create standard deviations NetCDF file.

s = c_std_bry(S);


%  Write out grid data.

v = 'spherical';    s = nc_write(S.ncname, v, S.spherical);
v = 'Vtransform';   s = nc_write(S.ncname, v, S.Vtransform);
v = 'Vstretching';  s = nc_write(S.ncname, v, S.Vstretching);

v = 'theta_s';      f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 'theta_b';      f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 'Tcline';       f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 'hc';           f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 's_rho';        f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 's_w';          f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 'Cs_r';         f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 'Cs_w';         f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);
v = 'h';            f = nc_read(INIfile,v);  s = nc_write(S.ncname,v,f);

v = 'lon_rho';      s = nc_write(S.ncname, v, S.rlon);
v = 'lat_rho';      s = nc_write(S.ncname, v, S.rlat);
v = 'lon_u';        s = nc_write(S.ncname, v, S.ulon);
v = 'lat_u';        s = nc_write(S.ncname, v, S.ulat);
v = 'lon_v';        s = nc_write(S.ncname, v, S.vlon);
v = 'lat_v';        s = nc_write(S.ncname, v, S.vlat);

if (S.curvilinear),
  v='angle';        s = nc_write(S.ncname, v, S.angle);
end,
  
if (S.masking),
  v='mask_rho';     s = nc_write(S.ncname, v, S.rmask);
  v='mask_u';       s = nc_write(S.ncname, v, S.umask);
  v='mask_v';       s = nc_write(S.ncname, v, S.vmask);
end,
    
% Write out standard deviation data.

rec = 1;
  
s = nc_write(S.ncname, 'ocean_time', 0, rec);

s = nc_write(S.ncname, 'zeta_obc', S.zeta_std, rec);
s = nc_write(S.ncname, 'ubar_obc', S.ubar_std, rec);
s = nc_write(S.ncname, 'vbar_obc', S.vbar_std, rec);
s = nc_write(S.ncname, 'u_obc'   , S.u_std   , rec);
s = nc_write(S.ncname, 'v_obc'   , S.v_std   , rec);
s = nc_write(S.ncname, 'temp_obc', S.temp_std, rec);
s = nc_write(S.ncname, 'salt_obc', S.salt_std, rec);
