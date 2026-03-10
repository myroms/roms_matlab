function add_jerlov(ncfile, wtype)

%
% ADD_JERLOV:  Adds Jerlov water type index to a ROMS Grid NetCDF file
%
% add_Jerlov(ncfile, wtype)
%
% This function adds spatially-varying Jerlov water type to an existng
% ROMS Grid NetCDF file. The wtype_grid variable is used when the CPP
% option WTYPE_GRID is activated.
%
% Currently, the water type classification is based on Jerlov water
% type using a double exponential function for light absorption:
%
%    Array
%    Index   WaterType   Examples
%    -----   ---------   --------
%
%      1         I       Open Pacific
%      2         IA      Eastern Mediterranean, Indian Ocean
%      3         IB      Western Mediterranean, Open Atlantic
%      4         II      Coastal waters, Azores
%      5         III     Coastal waters, North Sea
%      6         1       Skagerrak Strait
%      7         3       Baltic
%      8         5       Black Sea
%      9         7       Dark coastal water
%
% On Input:
%
%    ncfile        GRID NetCDF file name (string)
%
%    wtype         Jerlov water type index (floating-point)
%

% git $Id$
%=========================================================================%
%  Copyright (c) 2002-2026 The ROMS Group                                 %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.md                            Hernan G. Arango      %
%=========================================================================%

%--------------------------------------------------------------------------
% Inquire grid NetCDF file.
%--------------------------------------------------------------------------

I = nc_inq(ncfile);

got.wtype_grid  = any(strcmp({I.Variables.Name}, 'wtype_grid'));

if (any(strcmp({I.Variables.Name}, 'spherical')))
  spherical = nc_read(ncfile, 'spherical');
  if (ischar(spherical))
    if (spherical == 'T' || spherical == 't')
      spherical = true;
    else
      spherical = false;
    end
  end
else
  spherical = true;
end

%--------------------------------------------------------------------------
%  If appropriate, define Jerlov water type index.
%--------------------------------------------------------------------------

append_vars = false;

ic= 0;

if (~got.wtype_grid)
  ic = ic + 1;
  S.Dimensions = I.Dimensions;
  S.Variables(ic) = roms_metadata('wtype_grid', spherical);
  append_vars = true;
else
  disp(['Variable "wtype_grid" already exists. Updating its value']);
end

if (append_vars)
  check_metadata(S);
  nc_append(ncfile, S);
end

%--------------------------------------------------------------------------
%  Write out sponge variables into GRID NetCDF file.
%--------------------------------------------------------------------------

Lr = I.Dimensions(strcmp({I.Dimensions.Name},'xi_rho' )).Length;
Mr = I.Dimensions(strcmp({I.Dimensions.Name},'eta_rho')).Length;

[Im,Jm] = size(wtype);
if (Im == Lr && Jm == Mr)
  nc_write(ncfile, 'wtype_grid', wtype);
else
  error([' ADD_JERLOV: size(wgrid_type) is different to Lr = ',         ...
         num2str(Lr), blanks(3), 'Mr = ', num2str(Mr)]);
end

return
