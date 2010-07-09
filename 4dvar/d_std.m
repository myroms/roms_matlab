%
%  D_STD:  Driver script to compute and write 4D-Var standard deviations.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  standard deviations which are used to convert modeled error correlations
%  to error covariances.
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
%  month will be appended as *i_jan.nc:

STDfile = 'wc13_std_i.nc';

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

S.grd_file = GRDfile;

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
  Nvar = 0;
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
    disp([' D_STD: error while processing, rec = ', num2str(rec)]);
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

%  Accumulate monthly variance fields.

  disp(' ');
  disp([ 'Computing standard deviation, month = ', num2str(m), ' ...']);
  disp(' ');

  for n=1:nfiles,

    ncfile = fullfile(HISdir, HISfile(n).name);

    time = nc_read(ncfile,'ocean_time');
    Nrec = length(time);

    for rec=1:Nrec,

      [year,month,day]=datevec(datenum(1968,5,23) + time(rec)/86400);
      
      if (month == m),

        mydate=datestr(datenum(1968,5,23) + time(rec)/86400);
        disp([ '*** Processing Variance: ', mydate]);

        try,
          S.zeta_std = S.zeta_std + ...
	               (nc_read(ncfile,'zeta',rec) - S.zeta_avg) .^ 2;
          S.ubar_std = S.ubar_std + ...
	               (nc_read(ncfile,'ubar',rec) - S.ubar_avg) .^ 2;
          S.vbar_std = S.vbar_std + ...
	               (nc_read(ncfile,'vbar',rec) - S.vbar_avg) .^ 2;
          S.u_std    = S.u_std    + ...
	               (nc_read(ncfile,'u'   ,rec) - S.u_avg   ) .^ 2;
          S.v_std    = S.v_std    + ...
	               (nc_read(ncfile,'v'   ,rec) - S.v_avg   ) .^ 2;
          S.temp_std = S.temp_std + ...
	               (nc_read(ncfile,'temp',rec) - S.temp_avg) .^ 2;
          S.salt_std = S.salt_std + ...
	               (nc_read(ncfile,'salt',rec) - S.salt_avg) .^ 2;

          Nvar = Nvar + 1;
        catch,
          disp([' D_STD: error while processing, rec = ', num2str(rec)]);
          disp(['        in file: ', ncfile]);
          return        
	end,
      
      end,
    
    end,
  
  end,

%  Compute standard deviations.

  S.zeta_std = sqrt(S.zeta_std ./ (Nvar - 1));
  S.ubar_std = sqrt(S.ubar_std ./ (Nvar - 1));
  S.vbar_std = sqrt(S.vbar_std ./ (Nvar - 1));
  S.u_std    = sqrt(S.u_std    ./ (Nvar - 1));
  S.v_std    = sqrt(S.v_std    ./ (Nvar - 1));
  S.temp_std = sqrt(S.temp_std ./ (Nvar - 1));
  S.salt_std = sqrt(S.salt_std ./ (Nvar - 1));

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
  
  s = nc_write(S.ncname, 'ocean_time', 0, rec);

  s = nc_write(S.ncname, 'zeta', S.zeta_std, rec);
  s = nc_write(S.ncname, 'ubar', S.ubar_std, rec);
  s = nc_write(S.ncname, 'vbar', S.vbar_std, rec);
  s = nc_write(S.ncname, 'u'   , S.u_std   , rec);
  s = nc_write(S.ncname, 'v'   , S.v_std   , rec);
  s = nc_write(S.ncname, 'temp', S.temp_std, rec);
  s = nc_write(S.ncname, 'salt', S.salt_std, rec);

% Process next month.

end,
