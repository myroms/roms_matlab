%
%  D_STD_UNBALANCED:  Driver script to compute and write 4D-Var unbalaced
%                     standard deviations.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  standard deviations when the balance operator is activated. They are
%  used to convert modeled error correlations to error covariances.
%
%  The first step is to run the model application for a period that is
%  long enough to compute meaningful circulation statistics, like mean
%  and standard deviations for all prognostic state variables: zeta, u,
%  v, temp, and salt.
%
%  As an example, we compute the 4D-Var standard deviations for the WC13
%  application from daily history files from 1/1/2000 to 12/25/2004.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%  Set standard deviation NetCDF file. The file name is edited and the
%  month will be appended as *iu_jan.nc:

STDfile = 'wc13_std_iu.nc';

mstr = {'jan', 'feb', 'mar', 'apr', 'may', 'jun', ...
        'jul', 'aug', 'sep', 'oct', 'nov', 'dec'};

%  Set input grif file.

my_root = '/home/arango/ocean/toms/repository/test';

GRDfile = fullfile(my_root, 'WC13/Data', 'wc13_grd.nc');

%  Set input history files (string cell structure).

my_root = '/home/arango/ocean/toms/repository/test';

HISdir  = fullfile(my_root, 'WC13/STD/Data');

HISfile = dir(fullfile(HISdir, 'wc*.nc'));

nfiles = length(HISfile);

%---------------------------------------------------------------------------
%  Compute monthly averages and standard deviations.
%---------------------------------------------------------------------------

%  Extract first history file name.

HisFile1 = fullfile(HISdir, HISfile(1).name);

%  Set NetCDF file creation parameters.

[vnames,nvars] = nc_vname(HisFile1);

S.curvilinear = false;
S.masking     = false;
S.Vtransform  = 1;
S.Vstretching = 1;

for n=1:nvars,
  name=deblank(vnames(n,:));
  switch name
    case 'spherical'
      S.spherical = nc_read(HisFile1, 'spherical');
      if (ischar(S.spherical)),
        if (S.spherical == 'T' | S.spherical == 't');
          S.spherical = true;
        else,
          S.spherical = false;
        end,
      end,
    case 'Vtransform'
      S.Vtransform  = nc_read(HisFile1, 'Vtransform');
    case 'Vstretching'
      S.Vstretching = nc_read(HisFile1, 'Vstretching');
    case 'angle'
      S.angle = nc_read(HisFile1, 'angle');
      S.curvilinear = true;
    case 'mask_rho'
      S.rmask = nc_read(HisFile1, 'mask_rho');
      S.umask = nc_read(HisFile1, 'mask_u');
      S.vmask = nc_read(HisFile1, 'mask_v');
      S.masking = true;
  end,
end,

S.title = 'California Current System, 1/3 degree resolution (WC13)';

S.grid_file = GRDfile;

[Lr,Mr,Nr] = size(nc_read(HisFile1, 'temp', 1));

S.Lm = Lr - 2;                      % number of interior RHO x-points
S.Mm = Mr - 2;                      % number of interior RHO y-points
S.N  = Nr;                          % number of vertical RHO levels

S.do_zeta = true;                   % free-surface
S.do_ubar = true;                   % vertically integrated u-momentum
S.do_vbar = true;                   % vertically integrated v-momentum
S.do_u    = true;                   % u-momentum
S.do_v    = true;                   % v-momentum
S.do_temp = true;                   % temperature
S.do_salt = true;                   % salinity

%  Read in grid.

S.rlon = nc_read(HisFile1, 'lon_rho');
S.rlat = nc_read(HisFile1, 'lat_rho');

S.ulon = nc_read(HisFile1, 'lon_u');
S.ulat = nc_read(HisFile1, 'lat_u');

S.vlon = nc_read(HisFile1, 'lon_v');
S.vlat = nc_read(HisFile1, 'lat_v');

%  Process month-by-month.

for m=1:12,

%  Set standard deviation file name.

  lstr = length(STDfile)-3;
  S.ncname = char(strcat(STDfile(1:lstr), '_', mstr(m), '.nc'));
  
%  Initialize mean and variance arrays.

  Navg = 0;
  rec  = 1;

  S.month = m;

  try,
    S.zeta_avg = zeros(size(nc_read(HisFile1, 'zeta', rec)));
    S.ubar_avg = zeros(size(nc_read(HisFile1, 'ubar', rec)));
    S.vbar_avg = zeros(size(nc_read(HisFile1, 'vbar', rec)));
    S.u_avg    = zeros(size(nc_read(HisFile1, 'u'   , rec)));
    S.v_avg    = zeros(size(nc_read(HisFile1, 'v'   , rec)));
    S.temp_avg = zeros(size(nc_read(HisFile1, 'temp', rec)));
    S.salt_avg = zeros(size(nc_read(HisFile1, 'salt', rec)));
  catch,
    disp([' D_STD: error while processing, rec = ',num2str(rec)]);
    didp(['        in file: ', HisFile1]);
    return
  end,

  S.zeta_std = S.zeta_avg;
  S.ubar_std = S.ubar_avg;
  S.vbar_std = S.vbar_avg;
  S.u_std    = S.u_avg;
  S.v_std    = S.v_avg;
  S.temp_std = S.temp_avg;
  S.salt_std = S.salt_avg;
 
  disp(' ');
  disp([ 'Computing mean fields, month = ', num2str(m), ' ...']);
  disp(' ');

%  Accumulate montly fields.

  for n=1:nfiles,

    ncfile = fullfile(HISdir, HISfile(n).name);

    time = nc_read(ncfile,'ocean_time');
    Nrec = length(time);

    for rec=1:Nrec,

      [year,month,day]=datevec(datenum(1968,5,23) + time(rec)/86400);
      
      if (month == m),

        mydate=datestr(datenum(1968,5,23) + time(rec)/86400);
        disp([ '*** Processing Averages: ', mydate]);

        try,
          S.zeta_avg = S.zeta_avg + nc_read(ncfile, 'zeta', rec);
          S.ubar_avg = S.ubar_avg + nc_read(ncfile, 'ubar', rec);
          S.vbar_avg = S.vbar_avg + nc_read(ncfile, 'vbar', rec);
          S.u_avg    = S.u_avg    + nc_read(ncfile, 'u'   , rec);
          S.v_avg    = S.v_avg    + nc_read(ncfile, 'v'   , rec);
          S.temp_avg = S.temp_avg + nc_read(ncfile, 'temp', rec);
          S.salt_avg = S.salt_avg + nc_read(ncfile, 'salt', rec);

          Navg = Navg + 1;
        catch,
          disp([' D_STD: error while processing, rec = ', num2str(rec)]);
          disp(['        in file: ', ncfile]);
          return
        end,
      
      end,
    
    end,
  
  end,

%  Compute monthly mean fields.

  S.zeta_avg = S.zeta_avg ./ Navg;
  S.ubar_avg = S.ubar_avg ./ Navg;
  S.vbar_avg = S.vbar_avg ./ Navg;
  S.u_avg    = S.u_avg    ./ Navg;
  S.v_avg    = S.v_avg    ./ Navg;
  S.temp_avg = S.temp_avg ./ Navg;
  S.salt_avg = S.salt_avg ./ Navg;


%---------------------------------------------------------------------------
%  Compute monthly unbalanced error covariance standard deviations.
%---------------------------------------------------------------------------

  disp(' ');
  disp(['   Computing unbalanced standard deviation, month = ', ...
        num2str(m), ' ...']);
  disp(' ');

%  Initialize balance operator structure, A.  See 'balance_driver.m' for
%  discussion about some of these parameters. The default values for the
%  parameters are set in 'ini_balance.m', the user can overwrite such
%  values here.

  if (exist('A')),
    clear('A');         % clear structure A
  end,
    
  A.Gname = GRDfile;    % application's grid

%  Set switch for the computation of balanced, baroclinic free-surface
%  using an elliptic equation (A.elliptic=true) or by integrating the
%  hydrostatic equation (A.elliptic=false). If false, the integration
%  is from bottom to surface (A.LNM_depth=0) or from level of no
%  motion to the surface (A.LNM_depth>0).

% A.elliptic = false;   % integrate hydrostatic equation
  A.elliptic = true;    % solve SSH elliptic equation

  if (A.elliptic),
%   A.Niter = 300;      % Number of elliptical solver iterations
  end,                  % (default 200)
  
  if (~A.elliptic),
    A.LNM_depth = 0;    % integrate from bottom to surface
%   A.LNM_depth = 500;  % integrate from z=-500 to surface
  end,   
    
%  Internal parameters for the computation of the balanced salinity
%  in terms of the temperature. These parameters are set in
%  "ini_balance.m"

 
% A.ml_depth=150;       % mixed-layer depth (m, positive)
%                       % (default 100)

% A.dTdz-min=0.0001;    % minimum dT/dz allowed (Celsius/m)
%                       % (default 0.001)

%  Process monthly fields.

  Nvar = 0;             % initialize variance counter

  for n=1:nfiles,

    ncfile = fullfile(HISdir, HISfile(n).name);

    time = nc_read(ncfile,'ocean_time');
    Nrec = length(time);

    for rec=1:Nrec,

      [year,month,day]=datevec(datenum(1968,5,23)+time(rec)/86400);
      
      if (month == m),

        mydate=datestr(datenum(1968,5,23) + time(rec)/86400);
        disp([ '*** Processing Variance: ', mydate]);

%  Set time record to use in the computation of the thermal expansion
%  and saline contraction coefficients in "ini_balance.m", which are
%  used to compute the balanced salinity in "s_balance.m".

        A.Hname = ncfile;               % history file to process

        A.HisTimeRec = rec;             % needed in balance_4dvar

%  Get basic state to use in the balance operator. Read in selected record
%  from History NetCDF file, "A.Hname". Only temperature and salinity are
%  needed in "balance_4dbar.m".

        try,
          A.zeta = nc_read(A.Hname, 'zeta', A.HisTimeRec);
          A.ubar = nc_read(A.Hname, 'vbar', A.HisTimeRec);
          A.vbar = nc_read(A.Hname, 'ubar', A.HisTimeRec);
          A.u    = nc_read(A.Hname, 'u'   , A.HisTimeRec);
          A.v    = nc_read(A.Hname, 'v'   , A.HisTimeRec);
          A.temp = nc_read(A.Hname, 'temp', A.HisTimeRec);
          A.salt = nc_read(A.Hname, 'salt', A.HisTimeRec);
        catch,
          disp([' D_STD: error while processing, rec = ', num2str(rec)]);
          disp(['        in file: ', ncfile]);
          return        
        end,

%  Compute state anomalies from computed time mean.  Notice that only the
%  "A.deltaT" is used to stablish the balance part of the other state
%  variables.

        A.deltaZ = A.zeta - S.zeta_avg;
        A.deltaU = A.u    - S.u_avg;
        A.deltaV = A.v    - S.v_avg;
        A.deltaT = A.temp - S.temp_avg;        % needed in balance_4dvar
        A.deltaS = A.salt - S.salt_avg;
 
%  Set first guess free-surface for elliptic equation.  It is only used
%  when A.elliptic = 1. Use basic state values.

        A.zeta_guess = A.zeta;                 % needed in balance_4dvar

%  Compute balanced state anomalies.

        [K]=balance_4dvar(A);
 
%  Compute unbalanced state. Notice that the balaced operator is K^(-1). 
%  The minus sign below accounts for the inverse operator.

        A.deltaZ_u = A.deltaZ - K.deltaZ_b;
        A.deltaU_u = A.deltaU - K.deltaU_b;
        A.deltaV_u = A.deltaV - K.deltaV_b;
        A.deltaS_u = A.deltaS - K.deltaS_b;

%  Accumulate unbalanced error covariance matrix variance.

        if (Nvar == 0),                         % initialize
          zeta_var = zeros(size(A.zeta));
          u_var    = zeros(size(A.u));
          v_var    = zeros(size(A.v));
          temp_var = zeros(size(A.temp));
          salt_var = zeros(size(A.salt));
        end,

        Nvar = Nvar + 1;

        zeta_var = zeta_var + A.deltaZ_u .^ 2;
        u_var    = u_var    + A.deltaU_u .^ 2;
        v_var    = v_var    + A.deltaV_u .^ 2;
        temp_var = temp_var + A.deltaT   .^ 2;
        salt_var = salt_var + A.deltaS_u .^ 2;
  
      end,

    end,

  end,

%  Compute unbalanced error covariance matrix standard deviations.

  A.zeta_std = sqrt(zeta_var ./ (Nvar - 1));
  A.u_std    = sqrt(u_var    ./ (Nvar - 1));
  A.v_std    = sqrt(v_var    ./ (Nvar - 1));
  A.temp_std = sqrt(temp_var ./ (Nvar - 1));
  A.salt_std = sqrt(salt_var ./ (Nvar - 1));

  A.ubar_std = zeros(size(A.ubar));       % Not used in 4D-Var
  A.vbar_std = zeros(size(A.vbar));       % 3D applications

%---------------------------------------------------------------------------
%  Write out standard deviation fields.
%---------------------------------------------------------------------------

  disp(' ');
  disp([ 'Writting out standard deviation, file = ', S.ncname]);
  disp(' ');

%  Create standard deviations NetCDF file.

  s = c_std(S);
  
%  Write out grid data.

  v = 'spherical';    s = nc_write(S.ncname, v, S.spherical);
  v = 'Vtransform';   s = nc_write(S.ncname, v, S.Vtransform);
  v = 'Vstretching';  s = nc_write(S.ncname, v, S.Vstretching);

  v = 'theta_s';      f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 'theta_b';      f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 'Tcline';       f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 'hc';           f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 's_rho';        f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 's_w';          f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 'Cs_r';         f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 'Cs_w';         f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);
  v = 'h';            f = nc_read(HisFile1,v);  s = nc_write(S.ncname,v,f);

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
  
  s = nc_write(S.ncname, 'ocean_time', 0 , rec);

  s = nc_write(S.ncname, 'zeta', A.zeta_std, rec);
  s = nc_write(S.ncname, 'ubar', A.ubar_std, rec);
  s = nc_write(S.ncname, 'vbar', A.vbar_std, rec);
  s = nc_write(S.ncname, 'u'   , A.u_std   , rec);
  s = nc_write(S.ncname, 'v'   , A.v_std   , rec);
  s = nc_write(S.ncname, 'temp', A.temp_std, rec);
  s = nc_write(S.ncname, 'salt', A.salt_std, rec);

% Process next month.

end,
