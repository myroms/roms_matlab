function nc_inq(fname)

%
% NC_INQ:  Inquire about the contents of a NetCDF file
%
% nc_inq(fname)
%
% This gets and prints the contents of a NetCDF file.  It displays the
% dimensions variables.
%
% On Input:
%
%    fname       NetCDF file name or URL file name (character string)
%
% Adapted from J.V. Mansbridge (CSIRO) "inqcdf.m" M-file.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Check the number of arguments.

if (nargin < 1),
  help nc_inq
  return
end

%  Use the appropriate interface it the input file is URL from an OpenDAP
%  server.

url=nc_url(fname);

if (url),
  nc_inq_java(fname);
else
  nc_inq_mexnc(fname);
end

return

function nc_inq_java(fname)

%
% NC_INQ_java:  Inquire about the contents of a NetCDF file
%
% nc_inq_java(fname)
%
% This gets and prints the contents of a NetCDF file.  It displays the
% dimensions variables. It uses SNCTOOLS function "nc_info".

%  Inquire information from URL NetCDF file.

Info=nc_info(fname); 

%  Report available dimensions.

ndims=length(Info.Dimension);

disp(' ')
disp(['Available dimensions and values:']);
disp(' ')

unlimited=0;

for i=1:ndims,
  dimnam=Info.Dimension(i).Name;
  dimsiz=Info.Dimension(i).Length;
  unlimited=unlimited+Info.Dimension(i).Unlimited;

  if (i > 9),
    s=[' '  int2str(i) ') ' dimnam ' = ' int2str(dimsiz)];
  else
    s=['  ' int2str(i) ') ' dimnam ' = ' int2str(dimsiz)];
  end,
  disp(s)
end,

if (unlimited == 0),
  disp(' ')
  disp ('     None of the dimensions is unlimited.')
else
  disp(' ')
  for i=1:ndims,
    if (Info.Dimension(i).Unlimited),
      dimnam=Info.Dimension(i).Name;
      s=['     ' dimnam ' is unlimited in length.'];
      disp(s)
    end    
  end
end

%  Report available variables.

nvars=length(Info.Dataset);

disp(' ')
disp(['Available Variables:']);
disp(' ')

for i=1:3:nvars

  stri=int2str(i);

  if (length(stri) == 1)
    stri=[ ' ' stri];
  end
  varnam=Info.Dataset(i).Name;
  s=[ '  ' stri ') ' varnam ];
  addit=26-length(s);
  for j=1:addit
    s=[ s ' '];
  end
   
  if (i < nvars)
    stri=int2str(i+2);
    if (length(stri) == 1)
      stri=[ ' ' stri];
    end
    varnam=Info.Dataset(i+1).Name;
    s=[ s '  ' stri ') ' varnam ];
    addit=52-length(s);
    for j=1:addit
      s=[ s ' '];
    end
  end 
   
  if (i < nvars - 1)
    stri=int2str(i+3);
    if (length(stri) == 1)
      stri=[ ' ' stri];
    end
    varnam=Info.Dataset(i+2).Name;
    s=[ s '  ' stri ') ' varnam ];
  end 
  disp(s)
end

return


function nc_inq_mexnc(fname)

%
% NC_INQ_MEXNC:  Inquire about the contents of a NetCDF file
%
% nc_inq_mexnc(fname)
%
% This gets and prints the contents of a NetCDF file.  It displays the
% dimensions variables. It uses MEXNC functions. Therefore, it cannot
% process a URL OpenDAP file.
%
  
% Open the NetCDF file.
  
[ncid,rcode]=mexnc('ncopen',fname,'nowrite');
if (rcode == -1),
  error([ 'NC_INQ: ncopen - unable to open file: ' fname])
end
disp('  ');
disp(['NetCDF file: ',fname]);

% Supress all error messages from NetCDF.

mexnc('setopts',0);

% Inquire about the contents of NetCDf file. Display information.

[ndims,nvars,ngatts,recdim,rcode]=mexnc('ncinquire',ncid);
if rcode == -1
  error([ 'mexnc: ncinquire - error while inquiring file: ' fname])
end

% Get and print out information about the dimensions.

disp(' ')
disp(['Available dimensions and values:']);
disp(' ')
for i=0:ndims-1,
  [dimnam,dimsiz,rcode]=mexnc('ncdiminq',ncid,i);
  if (i > 8),
    s=[' '  int2str(i+1) ') ' dimnam ' = ' int2str(dimsiz)];
  else
    s=['  ' int2str(i+1) ') ' dimnam ' = ' int2str(dimsiz)];
  end,
  disp(s)
end,

if (recdim <= 0),
  disp(' ')
  disp ('     None of the dimensions is unlimited.')
else
  [dimnam,dimsiz,rcode]=mexnc('ncdiminq',ncid,recdim);
  s=['     ' dimnam ' is unlimited in length.'];
  disp(' ')
  disp(s)
end,

% Get and print information about the variables.

disp(' ')
disp(['Available Variables:']);
disp(' ')

for i=0:3:nvars-1

  stri=int2str(i+1);

  if (length(stri) == 1)
    stri=[ ' ' stri];
  end
  [varnam,vartyp,nvdims,vdims,nvatts,rcode]=mexnc('ncvarinq',ncid,i);
  s=[ '  ' stri ') ' varnam ];
  addit=26-length(s);
  for j=1:addit
    s=[ s ' '];
  end
   
  if (i < nvars-1)
    stri=int2str(i+2);
    if (length(stri) == 1)
      stri=[ ' ' stri];
    end
    [varnam,vartyp,nvdims,vdims,nvatts,rcode]=mexnc('ncvarinq',ncid,i+1);
    s=[ s '  ' stri ') ' varnam ];
    addit=52-length(s);
    for j=1:addit
      s=[ s ' '];
    end
  end 
   
  if (i < nvars - 2)
    stri=int2str(i+3);
    if (length(stri) == 1)
      stri=[ ' ' stri];
    end
    [varnam,vartyp,nvdims,vdims,nvatts,rcode]=mexnc('ncvarinq',ncid,i+2);
    s=[ s '  ' stri ') ' varnam ];
  end 
  disp(s)
end

[rcode]=mexnc('ncclose', ncid);
if (rcode == -1),
  error(['** ERROR ** ncclose: rcode = ' num2str(rcode)])
end

return
