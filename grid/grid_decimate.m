function S = grid_decimate(Gfactor, Iname, Oname, varargin)

%
% GRID_DECIMATE:  Decimates a ROMS file into a coarser file
%
% S = grid_decimate(Gfactor, Iname, Oname, Lplot, Gout)
%
% Given a fine-resolution ROMS solution NetCDF file (Iname), this
% function creates a coarser-resolution file (Oname) by decimating
% the finer VARIABLES by the specified factor (Gfactor). To ensure
% that both fine and coarse grids at RHO points coincide at the
% domain boundary, we need:
%
%   The number fine RHO-grid points (0:L, 0:M) must be multiples
%   of Gfactor
%
%               rem(L, Gfactor) = 0
%               rem(M, Gfactor) = 0
%
% On Input:
%
%    Gfactor    Grid decimation factor (only 2,3,4 are supported)
%
%    Iname      Input  finer   NetCDF filename (string)
%
%    Oname      Output coarser NetCDF filename (string)
%
%    Lplot      Switch to plot grid decimation diagram (logical)
%
%    Gout       Coarser Grid NetCDF filename or structure (OPTIONAL)
%
% On Output:
%
%    S          Coaser resolution NetCDF structure
%

% svn $Id$
%=======================================================================%
%  Copyright (c) 2002-2025 The ROMS Group                               %
%    Licensed under a MIT/X style license                               %
%    See License_ROMS.md                            Hernan G. Arango    %
%=======================================================================%

% Initialize.

switch numel(varargin)
  case 0
    Lplot = false;
    Gout = [];
    got_grid = false;
  case 1
    Lplot = varargin{1};
    Gout = [];
    got_grid = false;
  case 2
    Lplot = varargin{1};
    Gout = varargin{2};
    if (~isstruct(Gout))
      G = get_roms_grid(Gout);
      Gname = Gout;
    else
      G = Gout;
      Gname = Gout.grid_name;
    end
    got_grid = true;
end

S = [];

% Check decimation factor.

legal = Gfactor > 1 && Gfactor < 5;
if (~any(legal))
  error([' DECIMATE: illegal decimation factor, Gfactor = ',          ...
         num2str(Gfactor)]);
end

% If applicable, get finer grid and variable structures.

F = get_roms_grid(Iname);

% Set data sampling indices (PSI-points).

[Lp,Mp] = size(F.h);

L = Lp-1;
M = Mp-1;

% Check if both number of grid point are even or odd.

Imultiple = rem(L, Gfactor) == 0;
Jmultiple = rem(M, Gfactor) == 0;

if (~Imultiple)
  error([' DECIMATE: The number grid-points in the X-direction, ',    ...
         ' L = ', num2str(L), ' must be multiple Gfactor = ',         ...
         num2str(Gfactor), ', rem(L, Gfactor) = ',                    ...
         num2str(rem(L,Gfactor))]);
end

if (~Jmultiple)
  error([' DECIMATE: The number grid-points in the X-direction, ',    ...
         ' L = ', num2str(M), ' must be multiple Gfactor = ',         ...
         num2str(Gfactor), ', rem(M, Gfactor) = ',                    ...
         num2str(rem(M,Gfactor))]);
end

% Get decimation indices at RHO-, U-, V-, and PSI-points.

Ir = unique([1:Gfactor:Lp Lp]);  Im = length(Ir);
Jr = unique([1:Gfactor:Mp Mp]);  Jm = length(Jr);

Iu = unique(1:Gfactor:L);
Ju = Jr;

Iv = Ir;
Jv = unique(1:Gfactor:M);

Ip = Iu;
Jp = Jv;

%----------------------------------------------------------------------
% Plot decimation grid diagrams.
%----------------------------------------------------------------------

if (Lplot)

% Fine and Coarse grids fractional coordinates.
  
  [YrF, XrF] = meshgrid(0.5:1:Mp-0.5, 0.5:1:Lp-0.5);      % RHO-points
  [YpF, XpF] = meshgrid(1.0:1:M     , 1.0:1:L     );      % PSI-points
  [YuF, XuF] = meshgrid(0.5:1:Mp-0.5, 1.0:1:L     );      % U-points
  [YvF, XvF] = meshgrid(1.0:1:M     , 0.5:1:Lp-0.5);      % V-points

  delta = 0.5 * Gfactor - 0.5;

  XrC = XrF(Ir, Jr);        YrC = YrF(Ir, Jr);
  XpC = XpF(Ip, Jp)+delta;  YpC = YpF(Ip, Jp)+delta;
  XuC = XuF(Iu, Ju)+delta;  YuC = YuF(Iu, Ju);
  XvC = XvF(Iv, Jv);        YvC = YvF(Iv, Jv)+delta;
  
% Fine grid physical boundary perimeter.
  
  Istr = 1;  Iend = L;
  Jstr = 1;  Jend = M;

  XboxF = [squeeze(XpF(Istr:Iend,Jstr));                            ...
          squeeze(XpF(Iend,Jstr+1:Jend))';                          ...
          squeeze(flipud(XpF(Istr:Iend-1,Jend)));                   ...
          squeeze(fliplr(XpF(Istr,Jstr:Jend-1)))'];

  YboxF = [squeeze(YpF(Istr:Iend,Jstr));                            ...
          squeeze(YpF(Iend,Jstr+1:Jend))';                          ...
          squeeze(flipud(YpF(Istr:Iend-1,Jend)));                   ...
          squeeze(fliplr(YpF(Istr,Jstr:Jend-1)))'];

% Coarse grid physical boundary perimeter.

  Istr = 1;  Iend = Im-1;
  Jstr = 1;  Jend = Jm-1;

  XboxC = [squeeze(XpC(Istr:Iend,Jstr));                            ...
          squeeze(XpC(Iend,Jstr+1:Jend))';                          ...
          squeeze(flipud(XpC(Istr:Iend-1,Jend)));                   ...
          squeeze(fliplr(XpC(Istr,Jstr:Jend-1)))'];

  YboxC = [squeeze(YpC(Istr:Iend,Jstr));                            ...
          squeeze(YpC(Iend,Jstr+1:Jend))';                          ...
          squeeze(flipud(YpC(Istr:Iend-1,Jend)));                   ...
          squeeze(fliplr(YpC(Istr,Jstr:Jend-1)))'];
  
% Grid makers colors: RHO, PSI, U, and V.

  Color = ['k', 'g', 'b', 'r'];

% Finer Grid, SouthWestern Corner.
  
  figure;
  
  gr = plot(XrF, YrF, 'k:', XrF', YrF', 'k:');
  hold on;

  gp = plot(XboxF, YboxF, 'k-', 'LineWidth', 2);
  
  fc(1) = plot(XrF(:), YrF(:), 's',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(1));
  fc(2) = plot(XpF(:), YpF(:), 'o',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(2));
  fc(3) = plot(XuF(:), YuF(:), 'd',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(3));
  fc(4) = plot(XvF(:), YvF(:), 'v',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(4));

  axis ([0 5 0 5])
  title('Finer Grid: SouthWestern Corner');
  xlabel(['Decimation Factor = ', num2str(Gfactor)]);
  legend(fc,'RHO', 'PSI', 'U', 'V', 'location', 'NorthEastOutside');
  hold off;
  print('finer_SWcorner.png', '-dpng', '-r300'); 

% Finer Grid, NorthEastern Corner.
  
  figure;
  
  gr = plot(XrF, YrF, 'k:', XrF', YrF', 'k:');
  hold on;

  gp = plot(XboxF, YboxF, 'k-', 'LineWidth', 2);
  
  fc(1) = plot(XrF(:), YrF(:), 's',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(1));
  fc(2) = plot(XpF(:), YpF(:), 'o',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(2));
  fc(3) = plot(XuF(:), YuF(:), 'd',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(3));
  fc(4) = plot(XvF(:), YvF(:), 'v',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(4));

  axis ([L-5 L+1 M-5 M+1])
  title('Finer Grid: NorthEastern Corner');
  xlabel(['Decimation Factor = ', num2str(Gfactor)]);
  legend(fc,'RHO', 'PSI', 'U', 'V', 'location', 'NorthEastOutside');
  hold off;
  print('finer_NEcorner.png', '-dpng', '-r300'); 
  
% Coarser Grid, SouthWestern Corner.
  
  figure;
  
  gr = plot(XrF, YrF, 'k:', XrF', YrF', 'k:');
  hold on;

  gu = plot(XuF, YuF, 'k:', XuF', YuF', 'k:');
  gv = plot(XvF, YvF, 'k:', XvF', YvF', 'k:');

  gf = plot(XboxF, YboxF, 'k-', 'LineWidth', 2);
  gc = plot(XboxC, YboxC, 'm-', 'LineWidth', 2);
  
  fc(1) = plot(XrC(:), YrC(:), 's',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(1));
  fc(2) = plot(XpC(:), YpC(:), 'o',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(2));
  fc(3) = plot(XuC(:), YuC(:), 'd',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(3));
  fc(4) = plot(XvC(:), YvC(:), 'v',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(4));

  axis ([0 Gfactor*3 0 Gfactor*3])
  title('Coarser Grid: SouthWestern Corner');
  xlabel(['Decimation Factor = ', num2str(Gfactor)]);
  legend(fc,'RHO', 'PSI', 'U', 'V', 'location', 'NorthEastOutside');
  hold off;
  print('coarser_SWcorner.png', '-dpng', '-r300'); 
  
% Coarser Grid, NorthEastern Corner.
  
  figure;
  
  gr = plot(XrF, YrF, 'k:', XrF', YrF', 'k:');
  hold on;

  gu = plot(XuF, YuF, 'k:', XuF', YuF', 'm:');
  gv = plot(XvF, YvF, 'k:', XvF', YvF', 'm:');

  gf = plot(XboxF, YboxF, 'k-', 'LineWidth', 2);
  gc = plot(XboxC, YboxC, 'm-', 'LineWidth', 2);
  
  fc(1) = plot(XrC(:), YrC(:), 's',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(1));
  fc(2) = plot(XpC(:), YpC(:), 'o',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(2));
  fc(3) = plot(XuC(:), YuC(:), 'd',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(3));
  fc(4) = plot(XvC(:), YvC(:), 'v',                                 ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', Color(4));

  axis ([L-Gfactor*3 L+1 M-Gfactor*3 M+1])
  title('Coarser Grid: NorthEastern Corner');
  xlabel(['Decimation Factor = ', num2str(Gfactor)]);
  legend(fc,'RHO', 'PSI', 'U', 'V', 'location', 'NorthEastOutside');
  hold off;
  print('coarser_NEcorner.png', '-dpng', '-r300'); 

end

%----------------------------------------------------------------------
% Create output NetCDF from strucuture.
%----------------------------------------------------------------------

% Create a coarser grid NetCDF structure from finer grid. Then, edit
% the dimensions.

C = nc_inq(Iname);

dnames = {C.Dimensions.Name};

if (any(contains(dnames, 'xi_rho')))
  C.Dimensions(strcmp({C.Dimensions.Name},'xi_rho' )).Length = length(Ir);
  C.Dimensions(strcmp({C.Dimensions.Name},'eta_rho')).Length = length(Jr);
end

if (any(contains(dnames, 'xi_psi')))
  C.Dimensions(strcmp({C.Dimensions.Name},'xi_psi' )).Length = length(Ip);
  C.Dimensions(strcmp({C.Dimensions.Name},'eta_psi')).Length = length(Jp);
end
  
if (any(contains(dnames, 'xi_u')))
  C.Dimensions(strcmp({C.Dimensions.Name},'xi_u'   )).Length = length(Iu);
  C.Dimensions(strcmp({C.Dimensions.Name},'eta_u'  )).Length = length(Ju);
end

if (any(contains(dnames, 'xi_v')))
  C.Dimensions(strcmp({C.Dimensions.Name},'xi_v'   )).Length = length(Iv);
  C.Dimensions(strcmp({C.Dimensions.Name},'eta_v'  )).Length = length(Jv);
end
  
if (any(contains(dnames, 'IorJ')))
  IorJ = max(Im,Jm);
  obc  = C.Dimensions(strcmp({C.Dimensions.Name},'obc_adjust')).Length;
  C.Dimensions(strcmp({C.Dimensions.Name},'IorJ'   )).Length = IorJ;
end

if (any(contains(dnames, 'frc_adjust')))
  frc = C.Dimensions(strcmp({C.Dimensions.Name},'frc_adjust')).Length;
end

% Check number of vertical levels.

Index = strcmp({C.Dimensions.Name},'s_rho');
if (any(Index))
  N = C.Dimensions(strcmp({C.Dimensions.Name},'s_rho')).Length;
else
  N = 1;
end

% Check size of the record/unlimited dimension and initialize it to zero.

Index = [C.Dimensions.Unlimited];

if (any(Index))
  RecDim = C.Dimensions(Index).Name;
  Nrec = C.Dimensions(strcmp({C.Dimensions.Name}, RecDim)).Length;
  C.Dimensions(strcmp({C.Dimensions.Name}, RecDim)).Length = 0;
else
  Index = contains({C.Dimensions.Name}, 'time');
  RecDim = C.Dimensions(Index).Name;
  Nrec = C.Dimensions(strcmp({C.Dimensions.Name}, RecDim)).Length;
  C.Dimensions(strcmp({C.Dimensions.Name}, RecDim)).Length = 0;
end

% Report.

disp(blanks(1));
disp('Decimating ROMS Finer Grid NetCDF Data: ');
disp(blanks(1));
disp(['Input  File = ', Iname]);
disp(['Output File = ', Oname]);
if (got_grid)
  disp(['Output File = ', Gname]);
end  
disp(['Input  Grid = ', num2str(Lp),' x ',num2str(Mp),' x ',num2str(N)]);
disp(['Output Grid = ', num2str(Im),' x ',num2str(Jm),' x ',num2str(N)]);
disp(['    Gfactor = ', num2str(Gfactor)]);
disp(blanks(1));

% Create output NetCDF file.

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

nc_create(Oname, mode, C);

% Set global attributes.

status = nc_attadd(Oname, 'parent_file', Iname);
if (status ~= 0), return, end

status = nc_attadd(Oname, 'decimate_factor', int32(Gfactor));
if (status ~= 0), return, end

status = nc_attdel(Oname, 'history');
if (status ~= 0), return, end

history = ['File created using Matlab script: decimate, ',            ...
           date_stamp];
status = nc_attadd(Oname, 'history', history);
if (status ~= 0), return, end

%------------------------------------------------------------------------
% Process field data.
%------------------------------------------------------------------------

GridVars = {'h', 'f', 'pm', 'pn', 'angle',                            ...
            'lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',               ...
            'lon_u', 'lat_u', 'lon_v', 'lat_v',                       ...
            'mask_rho', 'mask_psi', 'mask_u', 'mask_v'};

nvars = length(C.Variables);

% First, process time-independent variables.

for n = 1:nvars
  Vname  = C.Variables(n).Name;
  nvdims = length(C.Variables(n).Dimensions);
  if (nvdims > 0)
    Dnames = {C.Variables(n).Dimensions.Name};
    RecVar = any(contains(Dnames, RecDim));

    if (RecVar)                           % Skip to the next variable;
      continue                            % time-dependent fields
    end                                   % are processed separately

    if (got_grid && any(contains(GridVars, Vname)))
      continue                            % Skip to the next variable;
    end                                   % get metrics from coarser file

    disp(['Processing time-independent variable: ', Vname]);
    
    Vtype  = 'NONE';
    if (any(strcmp(Dnames,'xi_rho'))),     Vtype = 'RHO'; end
    if (any(strcmp(Dnames,'xi_psi'))),     Vtype = 'PSI'; end
    if (any(strcmp(Dnames,'xi_u'))),       Vtype = 'U';   end
    if (any(strcmp(Dnames,'xi_v'))),       Vtype = 'V';   end
    if (any(strcmp(Dnames,'obc_adjust'))), Vtype = 'OBC'; end
    if (any(strcmp(Dnames,'frc_adjust'))), Vtype = 'FRC'; end
    Km = 0;
    if (any(strcmp(Dnames,'s_rho'))), Km = N; end
    if (any(strcmp(Dnames,'s_w'  ))), Km = N+1; end
      
    F = nc_read(Iname, Vname);            % if got_grid = false
    switch Vtype
      case ('RHO')                        % RHO-points metrics
        Fout = F(Ir, Jr);
        nc_write(Oname, Vname, Fout);
      case ('PSI')                        % PSI-points metrics
        Fout = F(Ip, Jp);
        nc_write(Oname, Vname, Fout);
      case ('U')                          % U-points metrics
        Fout = F(Iu, Ju);
        nc_write(Oname, Vname, Fout);
      case ('V')                          % V-points metrics
        Fout = F(Iv, Jv);
        nc_write(Oname, Vname, Fout);
      otherwise                           % Information arrays
        nc_write(Oname, Vname, F);
    end      
  else
    disp(['Processing time-independent variable: ', Vname]);

    F = nc_read(Iname, Vname);            % Process scalars
    nc_write(Oname, Vname, F);
  end
end

% Second, if a coarser grid is provided, write metric variables.

if (got_grid)
  for var = GridVars
    field = char(var);
    if (any(strcmp({C.Variables.Name}, field)))
      disp(['Processing grid metrics variable: ', Vname]);
      nc_write(Oname, field, G.(field));
    end
  end
end

% Finally, process time-dependent variables with record dimension.

for rec=1:Nrec    

  for n = 1:nvars
    Vname  = C.Variables(n).Name;
    nvdims = length(C.Variables(n).Dimensions);

    if (nvdims > 0)
      Dnames = {C.Variables(n).Dimensions.Name};
      RecVar = any(contains(Dnames, RecDim));

      if (RecVar)                         % time-dependent variables

        Vtype  = 'NONE';
        if (any(strcmp(Dnames,'xi_rho'))),     Vtype = 'RHO'; end
        if (any(strcmp(Dnames,'xi_psi'))),     Vtype = 'PSI'; end
        if (any(strcmp(Dnames,'xi_u'))),       Vtype = 'U';   end
        if (any(strcmp(Dnames,'xi_v'))),       Vtype = 'V';   end
        if (any(strcmp(Dnames,'obc_adjust'))), Vtype = 'OBC'; end
        if (any(strcmp(Dnames,'frc_adjust'))), Vtype = 'FRC'; end
        Km = 0;
        if (any(strcmp(Dnames,'s_rho'))), Km = N; end
        if (any(strcmp(Dnames,'s_w'  ))), Km = N+1; end
      
        disp(['Processing time-dependent variable: ', Vname, ', (',   ...
              Vtype, '-points)']);

        F = nc_read(Iname, Vname, rec);    
        if (nvdims > 2)                   % 2D or 3D fields
          switch Vtype
            case ('RHO')
              Ic = Ir;
              Jc = Jr;
            case ('PSI')
              Ic = Ip;
              Jc = Jp; 
            case ('U')
              Ic = Iu;
              Jc = Ju; 
            case ('V')
              Ic = Iv;
              Jc = Jv;       
            case ('OBC')
              if (Im >= Jm)
                Bc = Ir;
              else
                Bc = Jr;
              end
          end

          if (strcmp(Vtype, 'FRC'))       % Process FRC adjustment
            Im = length(Ic);
            Jm = length(Jc);
            Fout = zeros([Im Jm frc]);
            for k = 1:frc
              Fout(:,:,k) = squeeze(F(Ic,Jc,k));
            end    
            nc_write(Oname, Vname, Fout, rec);
          elseif (strcmp(Vtype, 'OBC'))   % Process OBC adjustment
            if (Km > 0)
              Fout = zeros([IorJ Km 4 obc]);
             for ibry = 1:4               % HGA need efficient code
                for iobc = 1:obc
                   for k = 1:Km
                     Fout(1:IorJ,k,ibry,iobc) = squeeze(F(Bc,k,ibry,iobc));
                   end
                end
              end
            else
              Fout = zeros([IorJ 4 obc]);
              for ibry = 1:4
                for iobc = 1:obc
                  Fout(1:IorJ,ibry,iobc) = squeeze(F(1:IorJ,ibry,iobc));
                end
              end
            end
            nc_write(Oname, Vname, Fout, rec);
          else          
            if (Km > 0)                   % Process 3D state field
              Im = length(Ic);
              Jm = length(Jc);
              Fout = zeros([Im Jm Km]);
              for k = 1:Km
                Fout(:,:,k) = squeeze(F(Ic,Jc,k));
              end  
            else                          % Process 2D state field
              Fout = F(Ic, Jc);
            end
            nc_write(Oname, Vname, Fout, rec);
          end  
        else                              % Process time coordinate
          nc_write(Oname, Vname, F, rec);
        end
      end
    end
  end
end

%------------------------------------------------------------------------
% Get full extracted grid structure.
%------------------------------------------------------------------------

S = get_roms_grid(Oname);

return
