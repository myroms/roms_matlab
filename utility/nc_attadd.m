function [status]=nc_attadd(fname,aname,avalue,vname);

%
% NC_ATTADD:  Add/modify a global or variable NetCDF attribute
%
% [status]=nc_attadd(fname,aname,avalue,vname);
%
% This function adds/modify a global or variable in a NetCDF file.
% If the "vname" argument is missing, it is assumed that "aname"
% is a global attribute.
%
% On Input:
%
%    fname      NetCDF file name (character string)
%    aname      Attribute name (character string)
%    avalue     Attribute value (numeric or character string)
%    vname      Variable name (character string; optional)
%
% On Output:
%
%    status     error flag
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%  Initialize error flag.

status=-1;

if (nargin < 4),
  got_var=0;
else,
  got_var=1;
end,

latt=length(aname);

%  Set NetCDF parameters.

[ncdouble]=mexnc('parameter','nc_double');
[ncfloat ]=mexnc('parameter','nc_float');
[ncglobal]=mexnc('parameter','nc_global');
[ncchar  ]=mexnc('parameter','nc_char');
[ncint   ]=mexnc('parameter','nc_int');

%  Open NetCDF file.

[ncid]=mexnc('open',fname,'nc_write');
if (ncid < 0),
  disp('  ');
  error(['NC_ATTADD: open - unable to open file: ', fname]);
  return
end

%  Put open file into define mode.

  [status]=mexnc('redef',ncid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    [status]=mexnc('close',ncid);
    error(['NC_ATTADD: redef - unable to put in definition mode: ',fname]);
    return
  end,

%---------------------------------------------------------------------------
%  Add/modify a variable attribute.
%---------------------------------------------------------------------------

if (got_var),

%  Get variable ID.

  [varid]=mexnc('inq_varid',ncid,vname);
  if (varid < 0),
    [status]=mexnc('close',ncid);
    disp('  ');
    error(['NC_ATTADD: inq_varid - cannot find variable: ',vname]);
    return
  end,

%  Inquire number of variable attributes.

  [nvatts,status]=mexnc('inq_varnatts',ncid,varid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_ATTADD: inq_varnatts - unable to inquire number of variable ' ...
	   'attributes: ',vname]);
  end,

%  Inquire if variable attribute exist.

  found=0;
  
  for i=0:nvatts-1
    [attnam,status]=mexnc('inq_attname',ncid,varid,i);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD: inq_attname: error while inquiring attribute: ' ...
	     num2str(i)]);
    end,
    lstr=length(attnam);
    if (strmatch(aname(1:latt),attnam(1:lstr),'exact')),
      found=1;
      break
    end,
  end,
  if (found),
    disp('  ');
    disp(['Requested attribute "',aname,'" already exist in variable "', ...
          vname,'".']);
    if ischar(avalue),
      disp(['Updating its value to: "',avalue,'".']);
    else,
      disp(['Updating its value to: ',num2str(avalue),'.']);
    end,
  end,

%  Inquire variable type.

  [vtype,status]=mexnc('inq_vartype',ncid,varid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_ATTADD: inq_vartype - unable to inquire datatype for ' ...
	   'variable: ',vname]);
  end,


%  Add/modify variable attribute.

  if ischar(avalue),
    lstr=length(avalue);
    [status]=mexnc('put_att_text',ncid,varid,aname,ncchar,lstr,avalue);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD: put_att_text - unable to define attribute ',...
             '"',aname,'"in variable: ',vname,'.']);
    end,
  else, 
    nval=length(avalue);
    switch (vtype)
      case (ncint)
        value=int32(avalue);
        [status]=mexnc('put_att_int',   ncid,varid,aname,vtype,nval,value);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_ATTADD: put_att_int - unable to define attribute ',...
                 '"',aname,'"in variable: ',vname,'.']);
        end,
      case (ncfloat)
        value=single(avalue);
        [status]=mexnc('put_att_float', ncid,varid,aname,vtype,nval,value)
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_ATTADD: put_att_float - unable to define attribute ',...
                 '"',aname,'"in variable: ',vname,'.']);
        end,
      case (ncdouble)
        value=double(avalue);
        [status]=mexnc('put_att_double',ncid,varid,aname,vtype,nval,value);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_ATTADD: put_att_double - unable to define attribute ',...
                 '"',aname,'"in variable: ',vname,'.']);
        end,
    end,
  end,

%---------------------------------------------------------------------------
%  Add/modify a global attribute.
%---------------------------------------------------------------------------
   
else,

%  Inquire number of global attributes.

  [natts,status]=mexnc('inq_natts',ncid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_ATTADD: inq_natts - unable to inquire number of global ' ...
	   'attributes: ',fname]);
  end,
  
%  Inquire if requested global attribute exist.

  found=0;
  
  for i=0:natts-1
    [attnam,status]=mexnc('inq_attname',ncid,ncglobal,i);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD: inq_attname: error while inquiring attribute: ' ...
	     num2str(i)]);
      return
    end,
    lstr=length(attnam);
    if (strmatch(aname(1:latt),attnam(1:lstr),'exact')),
      found=1;
      break
    end,
  end,
  if (found),
    disp('  ');
    disp(['Requested attribute "',aname,'" already exist in file: ', ...
          fname,'.']);
    disp(['Updating its value to: "',avalue,'".']);
  end,
  
%  Add/modify character attribute.

  if ischar(avalue),
    lstr=length(avalue);
    [status]=mexnc('put_att_text',ncid,ncglobal,aname,ncchar,lstr,avalue);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_ATTADD: put_att_text - unable to define attribute ',...
               '"',aname,'"in file: ',fname,'.']);
    end,
  end,
  
end,

%  Exit definition mode.

[status]=mexnc('enddef',ncid);
if (status < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_ATTADD: enddef - unable to exit definition mode: ',fname]);
  return
end,

%  Close NetCDF file.

[cstatus]=mexnc('ncclose',ncid);
if (cstatus < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_ATTADD: ncclose - unable to close NetCDF file: ', fname]);
end

return
