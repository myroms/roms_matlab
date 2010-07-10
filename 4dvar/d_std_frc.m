%
%  D_STD_frc:  Driver script to compute and write 4D-Var surface forcing
%              standard deviations.
%
%  This a user modifiable script that can be used to prepare ROMS 4D-Var
%  surface forcing standard deviations which are used to convert modeled
%  error correlations to error covariances.
%
%  The first step is to run the model application with air/sea bulk
%  fluxes (BULK_FLUXES option) for a period that is long enough to
%  compute meaningful circulation statistics, like mean and standard
%  deviations for surface forcing state variables: sustr, svstr,
%  shflux, and ssflux.
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

STDfile = 'wc13_std_f.nc';

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
          S.spherical = 1;
        else,
          S.spherical = 0;
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

S.do_sustr  = true;                 % surface u-momentum stress
S.do_svstr  = true;                 % surface v-momentum stress
S.do_shflux = true;                 % surface net heat flux
S.do_ssflux = true;                 % surface salt flux (E-P)*SALT

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
    S.sustr_avg  = zeros(size(nc_read(HisFile1, 'sustr' ,rec)));
    S.svstr_avg  = zeros(size(nc_read(HisFile1, 'svstr' ,rec)));
    S.shflux_avg = zeros(size(nc_read(HisFile1, 'shflux',rec)));
    S.ssflux_avg = zeros(size(nc_read(HisFile1, 'ssflux',rec)));
  catch,
    disp([' D_STD: error while processing, rec = ', num2str(rec)]);
    didp(['        in file: ', HisFile1]);
    return
  end,

  S.sustr_std  = S.sustr_avg;
  S.svstr_std  = S.svstr_avg;
  S.shflux_std = S.shflux_avg;
  S.ssflux_std = S.ssflux_avg;

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
          S.sustr_avg  = S.sustr_avg  + nc_read(ncfile, 'sustr' , rec);
          S.svstr_avg  = S.svstr_avg  + nc_read(ncfile, 'svstr' , rec);
          S.shflux_avg = S.shflux_avg + nc_read(ncfile, 'shflux', rec);
          S.ssflux_avg = S.ssflux_avg + nc_read(ncfile, 'ssflux', rec);

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

  S.sustr_avg  = S.sustr_avg  ./ Navg;
  S.svstr_avg  = S.svstr_avg  ./ Navg;
  S.shflux_avg = S.shflux_avg ./ Navg;
  S.ssflux_avg = S.ssflux_avg ./ Navg;

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
          S.sustr_std  = S.sustr_std + ...
                         (nc_read(ncfile,'sustr' ,rec) - S.sustr_avg)  .^ 2;
          S.svstr_std  = S.svstr_std + ...
                         (nc_read(ncfile,'svstr' ,rec) - S.svstr_avg)  .^ 2;
          S.shflux_std = S.shflux_std + ...
                         (nc_read(ncfile,'shflux',rec) - S.shflux_avg) .^ 2;
          S.ssflux_std = S.ssflux_std + ...
                         (nc_read(ncfile,'ssflux',rec) - S.ssflux_avg) .^ 2;

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

  S.sustr_std  = sqrt(S.sustr_std  ./ (Nvar - 1));
  S.svstr_std  = sqrt(S.svstr_std  ./ (Nvar - 1));
  S.shflux_std = sqrt(S.shflux_std ./ (Nvar - 1));
  S.ssflux_std = sqrt(S.ssflux_std ./ (Nvar - 1));

%---------------------------------------------------------------------------
%  Write out surface forcing standard deviation fields.
%---------------------------------------------------------------------------

  disp(' ');
  disp([ 'Writting out standard deviation, file = ', S.ncname]);
  disp(' ');

%  Create surface standard deviations NetCDF file.

  s = c_std_frc(S);


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
    
% Write out surface forcing standard deviation data.

  rec = 1;
  
  s = nc_write(S.ncname, 'ocean_time', 0, rec);

  s = nc_write(S.ncname, 'sustr' , S.sustr_std , rec);
  s = nc_write(S.ncname, 'svstr' , S.svstr_std , rec);
  s = nc_write(S.ncname, 'shflux', S.shflux_std, rec);
  s = nc_write(S.ncname, 'ssflux', S.ssflux_std, rec);

% Process next month.

end,
