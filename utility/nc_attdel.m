function [status]=nc_attdel(fname,aname,vname);

%
% NC_ATTDEL:  Delete requested NetCDF attribute
%
% [status]=nc_attdel(fname,aname,vname);
%
% This function deletes requested global or variable in a NetCDF file.
% If the "vname" argument is missing, it is assumed that "aname" is a
% global attribute.
%
% On Input:
%
%    fname      NetCDF file name (character string)
%    aname      Attribute name (character string)
%    vname      Variable name (character string; optional)
%
% On Output:
%
%    status     error flag
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%  Initialize error flag.

status=-1;

if (nargin < 3),
  got_var=0;
else,
  got_var=1;
end,

latt=length(aname);

%  Open NetCDF file.

[ncid]=mexnc('open',fname,'nc_write');
if (ncid < 0),
  disp('  ');
  error(['NC_ATTDEL: open - unable to open file: ', fname]);
  return
end

%  Put open file into define mode.

  [status]=mexnc('redef',ncid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    [status]=mexnc('close',ncid);
    error(['NC_ATTDEL: redef - unable to put in definition mode: ',fname]);
    return
  end,

%---------------------------------------------------------------------------
%  Deleting a variable attribute.
%---------------------------------------------------------------------------

if (got_var),

%  Get variable ID.

  [varid]=mexnc('inq_varid',ncid,vname);
  if (varid < 0),
    [status]=mexnc('close',ncid);
    disp('  ');
    error(['NC_ATTDEL: inq_varid - cannot find variable: ',vname]);
    return
  end,

%  Inquire number of variable attributes.

  [nvatts,status]=mexnc('inq_varnatts',ncid,varid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_ATTDEL: inq_varnatts - unable to inquire number of variable ' ...
	   'attributes: ',vname]);
  end,

%  Delete requested variable attribute.

  found=0;
  
  for i=0:nvatts-1
   [attnam,status]=mexnc('inq_attname',ncid,varid,i);
   if (status < 0),
     disp('  ');
     disp(mexnc('strerror',status));
     error(['NC_ATTDEL: inq_attname: error while inquiring attribute: ' ...
	    num2str(i)]);
   end,
   lstr=length(attnam);
   if (strmatch(aname(1:latt),attnam(1:lstr),'exact')),
     status=mexnc('del_att',ncid,varid,attnam);
     if (status < 0)
       disp('  ');
       disp(mexnc('strerror',status));
       error(['NC_ATTDEL: del_att - error while deleting attribute: ' aname]);
     end,
     found=1;
     break
   end,
 end,
 if (~found),
   disp('  ');
   disp(['Requested attribute "',aname,'" not found in variable "', ...
         vname,'".']);
 end,
  
%---------------------------------------------------------------------------
%  Deleting a global attribute.
%---------------------------------------------------------------------------
   
else,

%  Inquire number of global attributes.

  [natts,status]=mexnc('inq_natts',ncid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_ATTDEL: inq_natts - unable to inquire number of global ' ...
	   'attributes: ',fname]);
  end,
  
%  Delete requested global attribute.

  found=0;
  [NCGLOBAL]=mexnc('parameter','nc_global');
  
  for i=0:natts-1
    [attnam,status]=mexnc('inq_attname',ncid,NCGLOBAL,i);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_ATTDEL: inq_attname: error while inquiring attribute: ' ...
	     num2str(i)]);
      return
    end,
    lstr=length(attnam);
    if (strmatch(aname(1:latt),attnam(1:lstr),'exact')),
      status=mexnc('del_att',ncid,NCGLOBAL,attnam);
      if (status < 0)
        disp('  ');
        disp(mexnc('strerror',status));
        error(['NC_ATTDEL: del_att - error while deleting attribute: ' aname]);
      end,
      found=1;
      break
    end,
  end,
  if (~found),
    disp('  ');
    disp(['Requested global attribute "',aname,'" not found in: ', fname]);
  end,
  
end,

%  Exit definition mode.

[status]=mexnc('enddef',ncid);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_ATTDEL: enddef - unable to exit definition mode: ',fname]);
  return
end,

%  Close NetCDF file.

[cstatus]=mexnc('ncclose',ncid);
if (cstatus < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_ATTDEL: ncclose - unable to close NetCDF file: ', fname]);
end

return
