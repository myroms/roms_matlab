function info = nc_varinfo(ncfile,varname,field)

%
% NC_VARINFO  Returns metadata about a specific NetCDF variable.
%
%   INFO = NC_VARINFO(NCFILE,VARNAME)
%
% Returns a metadata structure about the variable VARNAME in the netCDF
% file NCFILE.
%
%   INFO will have the following fields:
%
%     Name      - A string containing the name of the variable.
%     Datatype  - The datatype of the variable.
%     Unlimited - Either 1 if the variable has an unlimited dimension or 
%                   0 if not.
%     Dimension - a cell array with the names of the dimensions upon 
%                   which this variable depends.
%     Size      - Size of the variable.
%     Attribute - An array of structures corresponding to the attributes 
%                   defined for the specified variable.
%                         
%   INFO = NC_VARINFO(NCFILE,VARNAME,<'field'>)
%
% Returns only one of the above fields: Datatype, Unlimited, Dimension, 
% Size, Attribute. Handy for use in expressions.
%
% Each "Attribute" element is a struct itself and contains the following 
% fields:
%
%     Name      - A string containing the name of the attribute.
%     Datatype  - The datatype of the variable.
%     Value     - Value of the attribute.
%
% Adapted from SNCTOOLS:
%
%   https://mexcdf.svn.sourceforge.net/svnroot/mexcdf/snctools/trunk
%
%   This function was adapted from SNCTOOLS function "nc_getvarinfo" to
%   return the variable information with the original numerical native
%   precision.  The original function "nc_getvarinfo" unwisely converted
%   all the numerical values to double precision.  This is problematic
%   when dealing with the variable numerical attributes. Specially, when
%   using the NetCDF attributes:  "_FillValue", "missing_value",
%   "scale_factor", and "add_offset".  And to lesser extend when using
%   the attributes: "valid_min", "valid_max", and "valid_range".
%
%   This function it is self contained to avoid interactions with the
%   original SNCTOOLS functions. All function calls are private and
%   attached below.  Only the Java interface is used since it can also
%   process files on OpenDAP servers.  The original SNCTOOLS private
%   functions were renamed by removing the prefix "get".
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%
  
info = nc_varinfo_java(ncfile,varname);

if nargin > 2,
  info = info.(field); 
end

return

function Dataset = nc_varinfo_java(ncfile,varname)

% Java backend for NC_VARINFO.

import ucar.nc2.dods.*     
import ucar.nc2.*         
                           
close_it = true;

% Try it as a local file.  If not a local file, try as via HTTP,
% then as dods.

if isa(ncfile,'ucar.nc2.NetcdfFile')
  jncid = ncfile;
  close_it = false;
elseif isa(ncfile,'ucar.nc2.dods.DODSNetcdfFile')
  jncid = ncfile;
  close_it = false;
elseif exist(ncfile,'file')
  fid = fopen(ncfile);
  ncfile = fopen(fid);
  fclose(fid);
  jncid = NetcdfFile.open(ncfile);
else
  try 
    jncid = NetcdfFile.open ( ncfile );
    catch %#ok<CTCH>
      try
        jncid = snc_opendap_open(ncfile);
        catch %#ok<CTCH>
          error ( 'NC_VARINFO: fileOpenFailure ', ...
                  'Could not open ''%s'' with java backend.' , ncfile);
      end
  end
end

if isa(varname,'ucar.nc2.Variable')
  jvarid = varname;
else
  jvarid = jncid.findVariable(varname);
  if isempty(jvarid)
    close(jncid);
    error ( 'NC_VARINFO: badVariableName ', ...
            'Could not locate variable %s', varname );
  end
end

% All the details are hidden here because we need the exact same
% functionality in "nc_info".

Dataset = nc_varidinfo_java(jvarid);

% If we were passed a java file id, don't close it upon exit.

if close_it
  close ( jncid );
end

return

function Dataset = nc_varidinfo_java(jvarid)

% NC_VARIDINFO_JAVA:  returns metadata structure for a netcdf variable
%
% This function is private to snctools.  It is called by nc_info and
% nc_varinfo, and uses the java API.
%
% USAGE:   Dataset = nc_varidinfo_java(jvarid);
% 
% PARAMETERS:
%
% Input:
%     jvarid:  
%         of type ucar.nc2.dods.DODSVariable
% Output:
%     Dataset:
%         struct of variable metadata

Attribute = struct('Name','','Nctype',0,'Datatype','','Value',NaN);

Dataset = struct('Name','','Nctype','','Datatype','','Unlimited',false,...
                 'Dimension',{''},'Size',[],'Attribute',Attribute,...
                 'Chunking',[],'Shuffle',0,'Deflate',0);

Dataset.Name = char ( jvarid.getName() );

% Get the datatype, store as an integer.

datatype = char(jvarid.getDataType().toString());

switch ( datatype )
  case 'double'
    Dataset.Nctype = nc_double;
    Dataset.Datatype = datatype;
  case 'float'
    Dataset.Nctype = nc_float;
    Dataset.Datatype = 'single';
  case {'int','long'}
    Dataset.Nctype = nc_int;
    Dataset.Datatype = 'int32';
  case 'short'
    Dataset.Nctype = nc_short;
    Dataset.Datatype = 'int16';
  case 'String'
    Dataset.Nctype = 12;           % Apparently, DODSNetcdfFile returns
    Dataset.Datatype = 'string';   % 'String', while NetcdfFile returns 'char'
  case 'char'
    Dataset.Nctype = nc_char;
    Dataset.Datatype = 'char';
  case 'byte'
    Dataset.Nctype = nc_byte;
    Dataset.Datatype = 'int8';
  otherwise
    error ( 'NC_VARINFO: unhandledDatatype ', ...
            '%s:  unhandled datatype ''%s''\n', datatype );
end

% Determine if it is unlimited or not.

Dataset.Unlimited = double ( jvarid.isUnlimited() );

% Retrieve the dimensions.

dims = jvarid.getDimensions();
nvdims = dims.size();
Dimension = cell(1,nvdims);

for j = 1:nvdims
  theDim = jvarid.getDimension(j-1);
  Dimension{j} = char ( theDim.getName() );
end
Dataset.Dimension = Dimension;

% Get the size of the variable

if nvdims == 0
  Dataset.Size = 1;
else
  Size = jvarid.getShape();
  Dataset.Size = Size';
end

if nc_getpref('PRESERVE_FVD')
  Dataset.Dimension = fliplr(Dataset.Dimension);
  Dataset.Size = fliplr(Dataset.Size);
end

% Get the list of attributes.

j_att_list = jvarid.getAttributes();
Dataset.Attribute = nc_attsinfo_java(j_att_list);

return

function Attribute = nc_attsinfo_java(j_att_list)

% NC_ATTSINFO_JAVA:  returns metadata about netcdf attributes
%
% USAGE:  Attribute = nc_attsinfo_java(j_att_list);
%
% PARAMETERS:
%
% Input:
%     j_att_list:
%         Of type "java.util.ArrayList".  Each list member is of type
%         "ucar.nc2.Attribute"
% Output:
%     Attribute:
%         Structure array of attribute metadata.  The fields are 
%         
%         Name
%         Nctype (backwards compatibility)
%         Datatype
%         Value

j_att_iterator = j_att_list.listIterator();
j = 0;

Attribute = struct('Name','','Nctype',0,'Datatype','','Value',[]);

while 1
    
% This throws an exception when we've reached the end of the list.

  try
    jatt = j_att_iterator.next();
    catch %#ok<CTCH>
      break;
  end
    
  j = j + 1;
  Attribute(j) = nc_attinfo_java(jatt);

end

if j == 0
  Attribute = [];
end

return

function Attribute = nc_attinfo_java(jatt)

% NC_ATTINFO_JAVA:  return metadata about netcdf attribute

Attribute = struct('Name','','Nctype',0,'Datatype','','Value',[]);
    
Attribute.Name = char(jatt.getName());
    
datatype = char(jatt.getDataType().toString());

switch ( datatype )
  case 'double'
    Attribute.Nctype = 6; %#ok<*AGROW>
    Attribute.Datatype = 'double';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = double(values)';
        
  case 'float'
    Attribute.Nctype = 5;
    Attribute.Datatype = 'single';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = single(values)';
        
  case 'String'
    Attribute.Nctype = 12;
    Attribute.Datatype = 'string';
    shape = double(jatt.getLength);
    Attribute.Value = snc_pp_strings( jatt, jatt.getValues(), shape) ;
        
  case 'char'
    Attribute.Nctype = 2;
    Attribute.Datatype = 'char';
    Attribute.Value = char ( jatt.getStringValue());
        
  case 'byte'
    Attribute.Nctype = 1;
    Attribute.Datatype = 'int8';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int8(values)';
        
  case 'short'
    Attribute.Nctype = 3;
    Attribute.Datatype = 'int16';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int16(values)';
        
  case 'int'
    Attribute.Nctype = 4;
    Attribute.Datatype = 'int32';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int32(values)';
        
  case 'long'
    Attribute.Nctype = 4;
    Attribute.Datatype = 'int64';
        
    j_array = jatt.getValues();
    values = j_array.copyTo1DJavaArray();
    Attribute.Value = int64(values)';
        
  otherwise
    error('NC_VARINFO: unhandledDatatype  ', ...
          'Unhandled attribute datatype ''%s''\n', datatype );
end

return

function value = nc_getpref(name)

%   PREF = NC_GETPREF(NAME) returns the value of an SNCTOOLS preference.
%   
%   This routine should not be called directly.

persistent PRESERVE_FVD

if isempty(PRESERVE_FVD)
    PRESERVE_FVD = getpref('SNCTOOLS','PRESERVE_FVD',false);
end

if strcmp(name,'PRESERVE_FVD')
    value = PRESERVE_FVD;
else
    error('unrecognized input to NC_GETPREF');
end

function jncid = snc_opendap_open(ncfile)

% SNC_OPENDAP_OPEN Open a connection to an OPeNDAP URL.  If the URL is
% password-protected, we will have to coerce netcdf-java to supply the
% credentials.

import ucar.nc2.dods.*  

% Is it a username/password protected URL?

pat = '(?<protocol>https{0,1})://(?<username>[^:]+):(?<password>[^@]+)@(?<host>[^/]+)(?<relpath>.*)';
parts = regexp(ncfile,pat,'names');

if numel(parts) == 0
  jncid = DODSNetcdfFile(ncfile);
else           % SncCreds is a custom java class supplied with SNCTOOLS.
  credentials = SncCreds(parts.username,parts.password);
  client = ucar.nc2.util.net.HttpClientManager.init(credentials,'snctools');
  opendap.dap.DConnect2.setHttpClient(client);
  ucar.unidata.io.http.HTTPRandomAccessFile.setHttpClient(client);
  ucar.nc2.dataset.NetcdfDataset.setHttpClient(client);
    
  jncid = DODSNetcdfFile(ncfile);
end

return

function values = snc_pp_strings(jobj,jdata,shape)

% Post process NC_STRING data into cell arrays.

if isempty(jdata)
  values = {''};
  return;
elseif strcmp(version('-release'),'14')
  % In R14, we must use the 2.2.x release of java.  No access to the
  % "getObject" method.  Assuming a single-valued string.
  values = {char(jobj.getStringValue())};
  return;
end

% Java says that the variable is laid out in row-major order.

if numel(shape) == 1
  values = cell([1 shape]);
else
  values = cell(shape);
end

for j = 1:prod(shape)
  values{j} = jdata.getObject(j-1);
end

return

