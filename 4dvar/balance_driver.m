%
% This a template script showing how to compute the 4DVar balance operator.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

%---------------------------------------------------------------------------
% Set Aplication Grid, History, and Average NetCDF files.
%---------------------------------------------------------------------------
%
% The user has a lot of latitude here.  In "ini_balance.m", the Grid
% NetCDF is used to read in bathymetry (h), Coriolis parameter (f),
% curvilinear metrics (pm, pn), and land/sea masking arrays, if any.
% The history file is used to compute the depths (Zr and Zw) and to
% compute the thermal expansion and saline contraction coefficients
% which are used in "s_balance.m".  When computing the depths, a
% free-surface value of zero is used since we need the rest state
% depths.

 my_root = '/home/arango/ocean/toms/repository/Projects/WC13';

 A.Gname = fullfile(my_root, 'Data',    'wc13_grd.nc');
 A.Hname = fullfile(my_root, 'Forward', 'wc13_his.nc');
 A.Aname = fullfile(my_root, 'Forward', 'wc13_avg.nc');

% Set output error covariance standard deviation file.

 A.Sname = 'wc13_std_i.nc';
 
%---------------------------------------------------------------------------
% Set switch for computation of balance, baroclinic free-surface. You
% may want to over-write the number of elliptic solver iterations. The
% default value is set internally to 200 iterations in "ini_balance.m".
%---------------------------------------------------------------------------

%A.elliptic = 0;     % integrate hydrostatic equation
 A.elliptic = 1;     % solve SSH elliptic equation

%A.Niter = 300;      % Number of elliptical solver iterations

%---------------------------------------------------------------------------
% It is possible to over-write the following internal parameters for the
% computation of the balanced salinity in terms of the temperature. These
% parameters are set in "ini_balance.m"
%---------------------------------------------------------------------------

%A.ml_depth=150;     % mixed-layer depth (m, positive), [default 100]

%A.dTdz-min=0.0001;  % minimum dT/dz allowed (Celsius/m), [default 0.001]

%---------------------------------------------------------------------------
% Compute time averaged state.  Use Average NetCDF file, "A.Aname". This
% is the best file to compute the actual time-averaged state.
%
% The strategy here is to run a long simulation for your application to
% compute the seasonal or annual statistics.  Again, you have a lot of
% latitude here.  You may want to do monthly.
%---------------------------------------------------------------------------
%
% Get number of time records in average file. You have the choice to
% select the values of "Tstr" and "Tend", which are optional arguments
% to "average.m", to compute the desired time-average window.

 Nrec = length(nc_read(A.Aname,'ocean_time'));
 
 Tstr = 1;
 Tend = Nrec;
 
 A.zeta_avg = average(A.Aname, 'zeta', Tstr, Tend);
 A.u_avg    = average(A.Aname, 'u'   , Tstr, Tend);
 A.v_avg    = average(A.Aname, 'v'   , Tstr, Tend);
 A.temp_avg = average(A.Aname, 'temp', Tstr, Tend);
 A.salt_avg = average(A.Aname, 'salt', Tstr, Tend);

%===========================================================================
% Compute balanced error covariance standard deviation.
%===========================================================================

 disp(' ');
 disp(['   Computing unbalanced error covariance standard deviation ...']);
 disp(' ');

% Set time records to process from the History NetCDF file.

 Nrec = length(nc_read(A.Hname,'ocean_time'));

 Tstr = 1;
 Tend = Nrec;

 ic = 0;                                           % counter for variance
 
 for rec=Tstr:Tend,

% Set time record to use in the computation of the thermal expansion
% and saline contraction coefficients in "ini_balance.m", which are
% used to compute the balanced salinity in "s_balance.m".

   A.HisTimeRec = rec;                             % needed in balance_4dvar
   
%---------------------------------------------------------------------------
% Get basic state to use in the balance operator. Read in selected record
% from History NetCDF file, "A.Hname".
%---------------------------------------------------------------------------
%
% Set "FillValue" to zero.  Although this is now the default, it
% forces users to update "nc_read.m".

   FillValue=0;

% Only temperature and salinity are needed in "balance_4dbar.m".

   A.zeta = nc_read(A.Hname, 'zeta', A.HisTimeRec, FillValue);
   A.u    = nc_read(A.Hname, 'u'   , A.HisTimeRec, FillValue);
   A.v    = nc_read(A.Hname, 'v'   , A.HisTimeRec, FillValue);
   A.temp = nc_read(A.Hname, 'temp', A.HisTimeRec, FillValue);
   A.salt = nc_read(A.Hname, 'salt', A.HisTimeRec, FillValue);
 
%---------------------------------------------------------------------------
% Compute state anomalies from computed time mean.  Notice that only the
% "A.deltaT" is used to stablish the balance part of the other state
% variables.
%---------------------------------------------------------------------------

   A.deltaZ = A.zeta - A.zeta_avg;
   A.deltaU = A.u    - A.u_avg;
   A.deltaV = A.v    - A.v_avg;
   A.deltaT = A.temp - A.temp_avg;                 % needed in balance_4dvar
   A.deltaS = A.salt - A.salt_avg;
 
%---------------------------------------------------------------------------
% Set first guess free-surface for elliptic equation.  It is only used
% when A.elliptic = 1.
%---------------------------------------------------------------------------
%
% Use basic state values.

   A.zeta_guess = A.zeta;                          % needed in balance_4dvar

%---------------------------------------------------------------------------
% Compute balanced state anomalies.
%---------------------------------------------------------------------------

   [K]=balance_4dvar(A);
 
%---------------------------------------------------------------------------
% Compute unbalanced state.
%---------------------------------------------------------------------------
%
% Compute balaced operator, K^(-1).  The minus sign below accounts for
% the inverse operator.

   A.deltaZ_u = A.deltaZ - K.deltaZ_b;
   A.deltaU_u = A.deltaU - K.deltaU_b;
   A.deltaV_u = A.deltaV - K.deltaV_b;
   A.deltaS_u = A.deltaS - K.deltaS_b;

%---------------------------------------------------------------------------
% Accumulate unbalance error covariance matrix variance.
%---------------------------------------------------------------------------

   if (rec == Tstr),
     zeta_var = zeros(size(A.zeta));
     u_var    = zeros(size(A.u));
     v_var    = zeros(size(A.v));
     temp_var = zeros(size(A.temp));
     salt_var = zeros(size(A.salt));
   end,

   ic = ic + 1;
   zeta_var = zeta_var + A.deltaZ_u.^2;
   u_var    = u_var    + A.deltaU_u.^2;
   v_var    = v_var    + A.deltaV_u.^2;
   temp_var = temp_var + A.deltaT.^2;
   salt_var = salt_var + A.deltaS_u.^2;
  
 end,

%---------------------------------------------------------------------------
% Compute unbalance error covariance matrix standard deviations.
%---------------------------------------------------------------------------

 A.zeta_std = sqrt(zeta_var ./ (ic-1));
 A.u_std    = sqrt(u_var    ./ (ic-1));
 A.v_std    = sqrt(v_var    ./ (ic-1));
 A.temp_std = sqrt(temp_var ./ (ic-1));
 A.salt_std = sqrt(salt_var ./ (ic-1));

%---------------------------------------------------------------------------
% Write out unbalanced error covariance standard devaitions to a NetCDF
% file
%---------------------------------------------------------------------------

 disp(' ');
 disp(['   Writting unbalanced standard deviation to file: ',A.Sname]);
 disp(' ');

% Set standard deviation creation parameters.

 S.ncname=A.Sname;
 S.Lm=K.Lm;
 S.Mm=K.Mm;
 S.N=K.N;

 [vname,nvars]=nc_vname(A.Hname);

 for n=1:nvars
   name=deblank(vname(n,:));
   switch name
     case 'spherical'
       spherical=nc_read(A.Hname,'spherical');
       if (spherical == 'T'),
	 S.spherical=1;
       else,
	 S.spherical=0;
       end,
     case 'Vtransform'
       S.Vtransform=nc_read(A.Hname,'Vtransform');
     case 'Vstretching'
       S.Vstretching=nc_read(A.Hname,'Vstretching');
     case 'hc'
       S.hc=nc_read(A.Hname,'hc');
     case 'theta_s'
       S.theta_s=nc_read(A.Hname,'theta_s');
     case 'theta_b'
       S.theta_b=nc_read(A.Hname,'theta_b');
     case 'Tcline'
       S.Tcline=nc_read(A.Hname,'Tcline');
     case 'sc_r'
       S.s_rho=nc_read(A.Hname,'sc_r');
     case 's_rho'
       S.s_rho=nc_read(A.Hname,'s_rho');
     case 'sc_w'
       S.s_w=nc_read(A.Hname,'sc_w');
     case 's_w'
       S.s_w=nc_read(A.Hname,'s_w');
     case 'Cs_r'
       S.Cs_r=nc_read(A.Hname,'Cs_r');
     case 'Cs_w'
       S.Cs_w=nc_read(A.Hname,'Cs_w');
   
   end,
 end,

 if (~isfield(S,'Vtransform'));
   S.Vtransform=1;
 end,

 if (~isfield(S,'Vstretching'));
   S.Vstretching=1;
 end,
 
% Set switches of variables to create (0: no, 1: yes).

 S.do_zeta = 1;                % free-surface
 S.do_ubar = 0;                % 2D U-velocity
 S.do_vbar = 0;                % 2D V-velocity
 S.do_u    = 1;                % 3D U-velocity
 S.do_v    = 1;                % 3D V-velocity
 S.do_temp = 1;                % temperature
 S.do_salt = 1;                % salinity
 
% Create standard deviation NetCDF file.

 [status]=c_std(S);  

%  Write out grid variables.

 if (S.spherical),
   [status]=nc_write(A.Sname, 'spherical',   'T');
 else,
   [status]=nc_write(A.Sname, 'spherical',   'F');
 end,

 [status]=nc_write(A.Sname, 'Vtransform',  S.Vtransform);
 [status]=nc_write(A.Sname, 'Vstretching', S.Vstretching);
 [status]=nc_write(A.Sname, 'theta_s',     S.theta_s);
 [status]=nc_write(A.Sname, 'theta_b',     S.theta_b);
 [status]=nc_write(A.Sname, 'Tcline',      S.Tcline);
 [status]=nc_write(A.Sname, 'hc',          S.hc);

 [status]=nc_write(A.Sname, 's_rho',       S.s_rho);
 [status]=nc_write(A.Sname, 's_w',         S.s_w);
 [status]=nc_write(A.Sname, 'Cs_r',        S.Cs_r);
 [status]=nc_write(A.Sname, 'Cs_w',        S.Cs_w);

 var = 'h';         f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);

 if (S.spherical),
   var = 'lon_rho'; f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'lat_rho'; f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'lon_u';   f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'lat_u';   f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'lon_v';   f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'lat_v';   f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
 else,
   var = 'x_rho';   f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'y_rho';   f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'x_u';     f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'y_u';     f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'x_v';     f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
   var = 'y_v';     f=nc_read(A.Gname,var); [status]=nc_write(A.Sname,var,f);
 end,

% Finally, write out unbalanced error covariance standard deviation fields.

 stdRec = 1;                               % NetCDF time record

 time=nc_read(A.Hname,'ocean_time',Tstr);  % standard deviation time (s)

 [status]=nc_write(A.Sname, 'ocean_time', time/86400, stdRec);

 if (S.do_zeta);
   [status]=nc_write(A.Sname, 'zeta', A.zeta_std, stdRec);
 end,
 
 if (S.do_ubar),
   [status]=nc_write(A.Sname, 'ubar', A.ubar_std, stdRec);
 end,

 if (S.do_vbar),
   [status]=nc_write(A.Sname, 'vbar', A.vbar_std, stdRec);
 end,

 if (S.do_u),
   [status]=nc_write(A.Sname, 'u',    A.u_std,    stdRec);
 end,

 if (S.do_v),
   [status]=nc_write(A.Sname, 'v',    A.v_std,    stdRec);
 end,

 if (S.do_temp),
   [status]=nc_write(A.Sname, 'temp', A.temp_std, stdRec);
 end,

 if (S.do_salt),
   [status]=nc_write(A.Sname, 'salt', A.salt_std, stdRec);
 end,
 
 disp(' ');
 disp('Done.');
 disp(' ');
