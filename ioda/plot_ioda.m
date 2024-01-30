function S = plot_ioda (ncname, varargin)
  
% PLOT_IODA:  Plots field from an IODA/UFO files.
%
% plot_ioda(ncname, ptype, wrtPNG)
%
% On Input:
%
%    ncname        IODA/UFO NetCDF4 filename (string)
%
%    group         NetCDF group values to plot (string, OPTIONAL)
%                    'hofx'        ~> Model at observation locations, H(x)
%                    'hofx0'       ~> Initial H(x)
%                    'hofx1'       ~> Final H(x)
%                    'ObsValue'    ~> observation value
%                    'oman'        ~> observation minus analysis
%                    'ombg'        ~> observation minus background
%
%    wrtPNG        Switch to write out PNG file (true or false)
%                    (Optional, default: false)
%
% On Output:
%
%    S             IODA/UFO data (struc)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2024 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

% Inquire input NetCDF file.

I = ncinfo(ncname);

Ngroups = length(I.Groups);

% Set optional arguments.

switch numel(varargin)
  case 0
    group   = "ObsValue";                        % default
    got_grp = false;
    wrtPNG  = false;
  case 1
    group   = varargin{1}
    if (isempty(group))
      group   = "ObsValue";
      got_grp = false;
    else
      got_grp = true; 
    end
    wrtPNG  = false;
  case 2
    group   = varargin{1};
    if (isempty(group))
      group   = "ObsValue";
      got_grp = false;
    else
      got_grp = true; 
    end
    wrtPNG  = varargin{2};
end
  
% Check if requested group is available

if (got_grp)
  if (~any(strcmp({I.Groups.Name}, group)))
    disp(blanks(1));
    disp(['NetCDF4 Group: ', group, ' not found in: ', ncname]);
    disp(blanks(1));
    disp('Available Groups:');
    disp(blanks(1));

    for i=1:3:Ngroups
      stri = int2str(i);
      if (length(stri) == 1)
        stri=[ ' ' stri];
      end
      grpnam = I.Groups(i).Name;
      s = [ '  ' stri ') ' grpnam ];
      addit = 26 - length(s);
      for j=1:addit
        s = [ s ' '];
      end
      if (i < Ngroups)
        stri = int2str(i+2);
        if (length(stri) == 1)
          stri = [ ' ' stri];
        end
        grpnam = I.Groups(i+1).Name;
        s = [ s '  ' stri ') ' grpnam ];
        addit = 52-length(s);
        for j=1:addit
          s = [ s ' '];
        end
      end 
      if (i < Ngroups - 1)
        stri = int2str(i+3);
        if (length(stri) == 1)
          stri = [ ' ' stri];
        end
        grpnam = I.Groups(i+2).Name;
        s = [ s '  ' stri ') ' grpnam ];
      end 
      disp(s);
    end
    disp(blanks(1));
    error(['PLOT_IODA: cannot find NetCDF4 Group: ', group]);  
  
  end  
end

% Read input IODA/UFO NetCDF4 file.

S = ioda_read(ncname);
M = ncinfo(ncname, 'MetaData');

% Set shortname (used for saving PNG files).

if (strfind(ncname, 'TS'))
  shortname = {'S', 'T'};
elseif (strfind(ncname, 'temp'))
  shortname = {'temp'};
elseif (strfind(ncname, 'salt'))
  shortname = {'salt'};
elseif (strfind(ncname, 'sst'))
  shortname = {'sst'};
elseif (strfind(ncname, 'sss'))
  shortname = {'sss'};
elseif (strfind(ncname, 'adt'))
  shortname = {'adt'};
elseif (strfind(ncname, 'usur'))
  shortname = {'usur'};
elseif (strfind(ncname, 'vsur'))
  shortname = {'vsur'};
elseif (strfind(ncname, 'uv_sur'))
  shortname = {'u_sur', 'v_sur'};
end 
S.shortname = shortname;

% Set Date for time axis labels.

t = S.datenum+S.dateTime/86400;
time = cellstr(datestr(t, 'dd-mmm-yyyy HH'));
Date = unique(cellstr(datestr(t, 'dd-mmm-yyyy HH')));

%--------------------------------------------------------------------------
% Plot data.
%--------------------------------------------------------------------------

LenFile  = length(ncname);
StrIndex = strfind(ncname, '/');
if (~isempty(StrIndex))
  MyFile = ncname(max(StrIndex)+1:LenFile);      % remove directory path
else
  MyFile = ncname;
end

if (isfield(S, 'sequenceNumber'))
  Nlocs = S.sequenceNumber;
else
  Nlocs = 1:S.nlocs;
end

% Plot for each state variable in the data structure.

for n = 1:S.nvars

  Vname = S.iodaVarName{n};

%--------------------------------------------------------------------------
% Plot "ObsValue" Group.  
%--------------------------------------------------------------------------

  if ( (~got_grp && any(strcmp({I.Groups.Name}, 'ObsValue'))) ||        ...
       (strcmp(group, 'ObsValue')) )

    if (any(strcmp({M.Variables.Name}, 'depth')))
      figure;
      h1 = plot(S.ObsValue{n}, S.depth, 'r+');
      xlabel(Vname);
      ylabel('depth');
      title(['Group:  ObsValue,   File:  ', untexlabel(MyFile)]);
      grid on;
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_obsZ.png');
        print(png_file, '-dpng', '-r300')
      end
    else
      figure;
      h2 = plot3(S.longitude, S.latitude, S.ObsValue{n}, 'r+');
      xlabel('Longitude');
      ylabel('latitude');
      zlabel(Vname);
      grid on;
      title(['Group:  ObsValue,   File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_obsXY.png');
        print(png_file, '-dpng', '-r300')
      end
    end  
  
  end
    
%--------------------------------------------------------------------------
% Plot "hofx" Group.  
%--------------------------------------------------------------------------

  if ( (~got_grp && any(strcmp({I.Groups.Name}, 'hofx'))) ||            ...
       (strcmp(group, 'hofx')) )

    figure;
    h1 = plot(Nlocs, S.ObsValue{n}, 'r+',                               ...
              Nlocs, S.hofx{n}, 'bo');
    legend('OBS', 'HofX', 'Location', 'best');
    xlabel('Observation Number')
    ylabel(Vname);
    grid on;
    title(['Group:  hofx,   File:  ', untexlabel(MyFile)]);
    if (wrtPNG)
      png_file = strcat(shortname{n}, '_enum.png');
      print(png_file, '-dpng', '-r300')
    end
       
    figure;
    h2 = plot(t, S.ObsValue{n}, 'r+',                                   ...
              t, S.hofx{n}, 'bo');
    legend('OBS', 'HofX', 'Location', 'best');
    datetick('x', 'dd-mmm-yyyy');
    xtickangle(0);
    ylabel(Vname);
    grid on;
    title(['Group:  hofx,   File:  ', untexlabel(MyFile)]);
    if (wrtPNG)
      png_file = strcat(shortname{n}, '_date.png');
      print(png_file, '-dpng', '-r300')
    end
       
    figure;
    h3 = scatter(Nlocs, S.ObsValue{n} - S.hofx{n}, 'bo',                ...
                 'MarkerFaceColor', 'b');
    hold on;
    line(Nlocs, zeros(size(Nlocs)), 'Color', 'red',                     ...
         'LineStyle', '--', 'LineWidth', 2);
    xlabel('Observation Number')
    ylabel([Vname, ' Innovation Vector']);
    title(['d = OBS - HofX for File:  ', untexlabel(MyFile)]);
    if (wrtPNG)
      png_file = strcat(shortname{n}, '_innov.png');
      print(png_file, '-dpng', '-r300')
    end
       
    if (any(strcmp({M.Variables.Name}, 'depth')))
      figure;
      h4 = plot(S.ObsValue{n}, S.depth, 'r+',                           ...
                S.hofx{n}, S.depth, 'bo');
      legend('OBS', 'HofX', 'Location', 'best');
      xlabel(Vname);
      ylabel('depth');
      grid on;
      title(['HofX File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_depth.png');
        print(png_file, '-dpng', '-r300')
      end
    else
      figure;
      h4 = plot3(S.longitude, S.latitude, S.ObsValue{n}, 'r+',         ...
                 S.longitude, S.latitude, S.hofx{n}, 'bo');
      legend('OBS', 'HofX', 'Location', 'best');
      xlabel('Longitude');
      ylabel('latitude');
      zlabel(Vname); 
      grid on;
      title(['Group:  hofx,   File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_lonlat.png');
        print(png_file, '-dpng', '-r300')
      end
    end
  
  end

%--------------------------------------------------------------------------
%  Plot "hofx0" Group: initial HofX.
%--------------------------------------------------------------------------
  
  if ( (~got_grp && any(strcmp({I.Groups.Name}, 'hofx0'))) ||           ...
       (strcmp(group, 'hofx0')) )

    figure;
    h1 = plot(Nlocs, S.ObsValue{n}, 'r+',                               ...
              Nlocs, S.hofx0{n}, 'bo');
    legend('OBS', 'Initial HofX', 'Location', 'best');
    xlabel('Observation Number')
    ylabel(Vname);
    grid on;
    title(['Group:  hofx0,   File:  ', untexlabel(MyFile)]);
    if (wrtPNG)
      png_file = strcat(shortname{n}, '_hofx0.png');
      print(png_file, '-dpng', '-r300')
    end
  
    if (any(strcmp({M.Variables.Name}, 'depth')))
      figure;
      h4 = plot(S.ObsValue{n}, S.depth, 'r+',                           ...
                S.hofx1{n}, S.depth, 'bo');
      legend('OBS', 'Initial HofX', 'Location', 'best');
      xlabel(Vname);
      ylabel('depth');
      grid on;
      title(['Group:  hofx0,   File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_hofx0_depth.png');
        print(png_file, '-dpng', '-r300')
      end
    else
      figure;
      h4 = plot3(S.longitude, S.latitude, S.ObsValue{n}, 'r+',         ...
                 S.longitude, S.latitude, S.hofx1{n}, 'bo');
      legend('OBS', 'Initial HofX', 'Location', 'best');
      xlabel('Longitude');
      ylabel('latitude');
      zlabel(Vname); 
      grid on;
      title(['Group:  hofx0,   File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_hofx0_lonlat.png');
        print(png_file, '-dpng', '-r300')
      end
    end
  
  end

%--------------------------------------------------------------------------
%  Plot "hofx1" Group: Final HofX.
%--------------------------------------------------------------------------
  
  if ( (~got_grp && any(strcmp({I.Groups.Name}, 'hofx1'))) ||           ...
       (strcmp(group, 'hofx1')) )

    figure;
    h1 = plot(Nlocs, S.ObsValue{n}, 'r+',                               ...
              Nlocs, S.hofx1{n}, 'bo');
    legend('OBS', 'Final HofX', 'Location', 'best');
    xlabel('Observation Number')
    ylabel(Vname);
    grid on;
    title(['Group:  hofx1,   File:  ', untexlabel(MyFile)]);
    if (wrtPNG)
      png_file = strcat(shortname{n}, '_hofx1.png');
      print(png_file, '-dpng', '-r300')
    end
 
    if (any(strcmp({M.Variables.Name}, 'depth')))
      figure;
      h4 = plot(S.ObsValue{n}, S.depth, 'r+',                           ...
                S.hofx1{n}, S.depth, 'bo');
      legend('OBS', 'Final HofX', 'Location', 'best');
      xlabel(Vname);
      ylabel('depth');
      grid on;
      title(['Group:  hofx1,   File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_hofx1_depth.png');
        print(png_file, '-dpng', '-r300')
      end
    else
      figure;
      h4 = plot3(S.longitude, S.latitude, S.ObsValue{n}, 'r+',         ...
                 S.longitude, S.latitude, S.hofx1{n}, 'bo');
      legend('OBS', 'Final HofX', 'Location', 'best');
      xlabel('Longitude');
      ylabel('latitude');
      zlabel(Vname); 
      grid on;
      title(['Group:  hofx1,   File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_hofx1_lonlat.png');
        print(png_file, '-dpng', '-r300')
      end
    end
  
  end

%--------------------------------------------------------------------------
%  Plot observation minus analysis
%--------------------------------------------------------------------------
  
  if ( (~got_grp && any(strcmp({I.Groups.Name}, 'oman'))) ||            ...
       (strcmp(group, 'oman')) )

    figure;
    h5 = bar(Nlocs, S.oman{n}, 'grouped', 'FaceColor', [.5, 0, .5]);
    xlabel('Observation Number')
    ylabel(Vname);
    grid on;
    title(['oman = OBS - Analysis,   File:  ', untexlabel(MyFile)]);    
    if (wrtPNG)
      png_file = strcat(shortname{n}, '_oman.png');
      print(png_file, '-dpng', '-r300')
    end
  
  end

%--------------------------------------------------------------------------
%  Plot observation minus background.
%--------------------------------------------------------------------------
  
  if ( (~got_grp && any(strcmp({I.Groups.Name}, 'oman'))) ||            ...
       (strcmp(group, 'ombg')) )

    figure;
    h5 = bar(Nlocs, S.ombg{n}, 'grouped', 'FaceColor', [0, .5, .5]);
    xlabel('Observation Number')
    ylabel(Vname);
    grid on;
    title(['ombg = OBS - Background,   File:  ', untexlabel(MyFile)]);    
    if (wrtPNG)
      png_file = strcat(shortname{n}, '_ombg.png');
      print(png_file, '-dpng', '-r300')
    end
  
  end
  
% plot next state variable, if any.
  
end

return
