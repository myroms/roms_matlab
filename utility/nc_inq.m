function nc_inq(fname)

%
% NC_INQ:  Inquire about the contents of a NetCDF file
%
% []=nc_inq(fname)
%
% This gets and prints the contents of a NetCDF file.  It displays the
% dimensions variables.
%
% On Input:
%
%    fname       NetCDF file name (character string)
%
% Adapted from J.V. Mansbridge (CSIRO) "inqcdf.m" M-file.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Check the number of arguments.

if (nargin < 1),
  help nc_inq
  return
end

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
