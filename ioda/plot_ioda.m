function S = plot_ioda (ncname, varargin)
  
% PLOT_IODA:  Plots field from an IODA/UFO files.
%
% plot_ioda(ncname, ptype, wrtPNG)
%
% On Input:
%
%    ncname        IODA/UFO NetCDF4 filename (string)
%
%    ptype         Plot type (string, OPTIONAL)
%                    'hofx'
%                    'obs'
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
%  Copyright (c) 2002-2023 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set optional arguments.

switch numel(varargin)
  case 0
    if (any(strfind(ncname, 'hofx')))
      ptype = 'hofx';
    else
      ptype = 'obs';
    end
    wrtPNG = false;
  case 1
    ptype = varargin{1};
    wrtPNG = false;
  case 2
    ptype = varargin{1};
    wrtPNG = varargin{2};
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

% Plot data.

LenFile  = length(ncname);
StrIndex = strfind(ncname, '_0000.')-1;
if (~isempty(StrIndex))
  MyFile = strcat(ncname(1:StrIndex), ncname(LenFile-3:LenFile));
else
  MyFile = ncname;
end
  
for n = 1:S.nvars

  Vname = S.iodaVarName{n};
  Nlocs = 1:S.nlocs;

  switch ptype
   
   case {'hofx', 'HofX'}

     if (isfield(S, 'hofx'))

       figure;
       h1 = plot(Nlocs, S.ObsValue{n}, 'r+',                            ...
                 Nlocs, S.hofx{n}, 'bo');
       legend('OBS', 'HofX', 'Location', 'best');
       xlabel('Observation Number')
       ylabel(Vname);
       grid on;
       title(['HofX File:  ', untexlabel(MyFile)]);
       if (wrtPNG)
	 png_file = strcat(shortname{n}, '_enum.png');
         print(png_file, '-dpng', '-r300')
       end
       
       figure;
       h2 = plot(t, S.ObsValue{n}, 'r+',                                ...
                 t, S.hofx{n}, 'bo');
       legend('OBS', 'HofX', 'Location', 'best');
       datetick('x', 'dd-mmm-yyyy');
       xtickangle(0);
       ylabel(Vname);
       grid on;
       title(['HofX File:  ', untexlabel(MyFile)]);
       if (wrtPNG)
	 png_file = strcat(shortname{n}, '_date.png');
         print(png_file, '-dpng', '-r300')
       end
       
       figure;
       h3 = scatter(Nlocs, S.ObsValue{n} - S.hofx{n}, 'bo',             ...
                    'MarkerFaceColor', 'b');
       hold on;
       line(Nlocs, zeros(size(Nlocs)), 'Color', 'red',                  ...
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
         h4 = plot(S.ObsValue{n}, S.depth, 'r+',                        ...
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
         h4 = plot3(S.longitude, S.latitude, S.ObsValue{n}, 'r+',       ...
                    S.longitude, S.latitude, S.hofx{n}, 'bo');
         legend('OBS', 'HofX', 'Location', 'best');
         xlabel('Longitude');
         ylabel('latitude');
         zlabel(Vname);
         grid on;
         title(['HofX File:  ', untexlabel(MyFile)]);
         if (wrtPNG)
	   png_file = strcat(shortname{n}, '_lonlat.png');
           print(png_file, '-dpng', '-r300')
         end
	 
       end
     end
   
   case {'obs', 'OBS'}

    if (any(strcmp({M.Variables.Name}, 'depth')))

      figure;
      h1 = plot(S.ObsValue{n}, S.depth, 'r+');
      xlabel(Vname);
      ylabel('depth');
      title(['OBS File:  ', untexlabel(MyFile)]);
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
      title(['OBS File:  ', untexlabel(MyFile)]);
      if (wrtPNG)
        png_file = strcat(shortname{n}, '_obsXY.png');
        print(png_file, '-dpng', '-r300')
      end
    
    end

  end
end

return
