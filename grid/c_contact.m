function c_contact(ncname, spherical, Ngrids, Ndatum)

%
% C_CONTACT:  Creates ROMS nesting Contact Points NetCDF file
%
% c_contact(ncname, spherical, Ngrids, Ndatum)
%
% This function creates a new Contact Points NetCDF file for a ROMS
% nesting application. The Contact Points data is written elsewhere.
% If the Ndatum argument is not provided, the "datum" dimension is
% is "unlimited".
%
% On Input:
%
%    ncname      Ouptut Contact Points NetCDF file name (string)
%
%    spherical   Spherical grid switch
%
%    Ngrids      Number of nested grids
%
%    Ndatum      Total number of contact points (optional)
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Set NetCDF file dimensions.
  
Ncontact = (Ngrids-1)*2;
Nweights = 4;

if (nargin < 4),
  Ndatum = netcdf.getConstant('UNLIMITED');
end

% ROMS FillValue.

FillValue = 1.0d+37;

%--------------------------------------------------------------------------
% Create Contact Points NetCDF.
%--------------------------------------------------------------------------

% Create file.

mode = netcdf.getConstant('CLOBBER');
mode = bitor(mode,netcdf.getConstant('64BIT_OFFSET'));

ncid = netcdf.create(ncname,mode);

% Define dimensions.

Did.Ngrids   = netcdf.defDim(ncid, 'Ngrids'  ,Ngrids);
Did.Ncontact = netcdf.defDim(ncid, 'Ncontact',Ncontact);
Did.Nweights = netcdf.defDim(ncid, 'Nweights',Nweights);
Did.datum    = netcdf.defDim(ncid, 'datum'   ,Ndatum);

% Define global attributes.

varid =  netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid, varid, 'type',                                      ...
              'ROMS Nesting Contact Regions Data');

%--------------------------------------------------------------------------
% Define variables.
%--------------------------------------------------------------------------

% Spherical switch.

varid = netcdf.defVar(ncid, 'spherical', nc_int, []);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'grid type logical switch') ;
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'Cartesian spherical');

% Number of interior RHO-points.

varid = netcdf.defVar(ncid, 'Lm', nc_int, Did.Ngrids);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'number of interior RHO-points in the I-direction') ;

varid = netcdf.defVar(ncid, 'Mm', nc_int, Did.Ngrids);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'number of interior RHO-points in the J-direction') ;

% Information variables.

varid = netcdf.defVar(ncid, 'coincident', nc_int, Did.Ngrids);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'coincident donor and receiver grids logical switch') ;
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'false true');

varid = netcdf.defVar(ncid, 'composite', nc_int, Did.Ngrids);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'composite grid type logical switch') ;
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'false true');

varid = netcdf.defVar(ncid, 'mosaic', nc_int, Did.Ngrids);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'mosaic grid type logical switch') ;
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'false true');

varid = netcdf.defVar(ncid, 'refinement', nc_int, Did.Ngrids);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'refinement grid type logical switch') ;
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'false true');

varid = netcdf.defVar(ncid, 'refine_factor', nc_int, Did.Ngrids);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'refinement factor from donor grid') ;

% Contact Points vertical interpolation switch.

varid = netcdf.defVar(ncid, 'interpolate', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'vertical interpolation at contact points logical switch');
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'false true');

% Contact region donor and receiver grid.

varid = netcdf.defVar(ncid, 'donor_grid', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'data donor grid number');

varid = netcdf.defVar(ncid, 'receiver_grid', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'data receiver grid number');

% Contact points start and end indices in data vector.

varid = netcdf.defVar(ncid, 'NstrR', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'starting contact RHO-point index in data vector');

varid = netcdf.defVar(ncid, 'NendR', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'ending contact RHO-point index in data vector');

varid = netcdf.defVar(ncid, 'NstrU', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'starting contact U-point index in data vector');

varid = netcdf.defVar(ncid, 'NendU', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'ending contact U-point index in data vector');

varid = netcdf.defVar(ncid, 'NstrV', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'starting contact V-point index in data vector');

varid = netcdf.defVar(ncid, 'NendV', nc_int, Did.Ncontact);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'ending contact V-point index in data vector');

% Contact region for each contact point.

varid = netcdf.defVar(ncid, 'contact_region', nc_int, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'contact region number');

% Contact points on receiver grid physical boundary.

varid = netcdf.defVar(ncid, 'on_boundary', nc_int, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
	      'contact point on receiver grid physical boundary') ;
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [int32(0), int32(1), int32(2), int32(3), int32(4)]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'false western southern eastern northern');

% Donor grid cell indices.

varid = netcdf.defVar(ncid, 'Idg', nc_int, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'I-left index of donor cell containing contact point');

varid = netcdf.defVar(ncid, 'Jdg', nc_int, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'J-bottom index of donor cell containing contact point');

% Receiver grid indices.

varid = netcdf.defVar(ncid, 'Irg', nc_int, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'receiver grid I-index of contact point');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');

varid = netcdf.defVar(ncid, 'Jrg', nc_int, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'receiver grid J-index of contact point');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');

% Receiver grid contact point coordinates.
  
varid = netcdf.defVar(ncid, 'Xrg', nc_double, Did.datum);
if (spherical),
  netcdf.putAtt(ncid, varid, 'long_name',                               ...
                'longitude of receiver grid contact points');
  netcdf.putAtt(ncid, varid, 'units',                                   ...
                'degree_east');
  netcdf.putAtt(ncid, varid, 'standard_name',                           ...
                'longitude');
else
  netcdf.putAtt(ncid, varid, 'long_name',                               ...
                'X-location of receiver grid contact points');
  netcdf.putAtt(ncid, varid, 'units',                                   ...
                'meter');
end

varid = netcdf.defVar(ncid, 'Yrg', nc_double, Did.datum);
if (spherical),
  netcdf.putAtt(ncid, varid, 'long_name',                               ...
                'latitude of receiver grid contact points');
  netcdf.putAtt(ncid, varid, 'units',                                   ...
                'degree_north');
  netcdf.putAtt(ncid, varid, 'standard_name',                           ...
                'latitude');
else
  netcdf.putAtt(ncid, varid, 'long_name',                               ...
                'Y-location of receiver grid contact points');
  netcdf.putAtt(ncid, varid, 'units',                                   ...
                'meter');
end

% Horizontal interpolation weights.

varid = netcdf.defVar(ncid, 'Hweight', nc_double,                       ...
                      [Did.Nweights Did.datum]);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'horizontal interpolation weights');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');

% Grid bathymetry.

varid = netcdf.defVar(ncid, 'h', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'bathymetry at RHO-points');
netcdf.putAtt(ncid, varid, 'units',                                     ...
              'meter');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');
netcdf.putAtt(ncid, varid, '_FillValue',                                ...
              FillValue);

% Coriolis parameter.

varid = netcdf.defVar(ncid, 'f', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'Coriolis parameter at RHO-points');
netcdf.putAtt(ncid, varid, 'units',                                     ...
              'second-1');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');
netcdf.putAtt(ncid, varid, '_FillValue',                                ...
              FillValue);

% Curvilinear metrics.

varid = netcdf.defVar(ncid, 'pm', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'curvilinear coordinate metric in XI');
netcdf.putAtt(ncid, varid, 'units',                                     ...
              'meter-1');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
             'Xrg Yrg');
netcdf.putAtt(ncid, varid, '_FillValue',                                ...
              FillValue);

varid = netcdf.defVar(ncid, 'pn', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'curvilinear coordinate metric in ETA');
netcdf.putAtt(ncid, varid, 'units',                                     ...
              'meter-1');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');
netcdf.putAtt(ncid, varid, '_FillValue',                                ...
              FillValue);

varid = netcdf.defVar(ncid, 'dndx', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'XI-derivative of inverse metric factor pn');
netcdf.putAtt(ncid, varid, 'units',                                     ...
              'meter');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');
netcdf.putAtt(ncid, varid, '_FillValue',                                ...
              FillValue);

varid = netcdf.defVar(ncid, 'dmde', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'ETA-derivative of inverse metric factor pm');
netcdf.putAtt(ncid, varid, 'units',                                     ...
              'meter');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');
netcdf.putAtt(ncid, varid, '_FillValue',                                ...
              FillValue);

% Grid rotation angle.

varid = netcdf.defVar(ncid, 'angle', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'angle between XI-axis and EAST');
netcdf.putAtt(ncid, varid, 'units',                                     ...
              'radians');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');
netcdf.putAtt(ncid, varid, '_FillValue',                                ...
              FillValue);

% Land/Sea masking.

varid = netcdf.defVar(ncid, 'mask', nc_double, Did.datum);
netcdf.putAtt(ncid, varid, 'long_name',                                 ...
              'land-sea mask of contact points');
netcdf.putAtt(ncid, varid, 'flag_values',                               ...
              [0 1]);
netcdf.putAtt(ncid, varid, 'flag_meanings',                             ...
              'land water');
netcdf.putAtt(ncid, varid, 'coordinates',                               ...
              'Xrg Yrg');

%--------------------------------------------------------------------------
% Close file.
%--------------------------------------------------------------------------

netcdf.endDef(ncid);

netcdf.close(ncid);

return
