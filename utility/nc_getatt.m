function [avalue]=nc_getatt(fname,aname,vname);

%
% NC_GETATT:  Gets a global or variable NetCDF attribute
%
% [avalue]=nc_getatt(fname,aname,vname)
%
% This function reads the spcedified global or variable attribute
% from input NetCDF file. If the "vname" argument is missing, it
% is assumed that "aname" is a global attribute.
%
% On Input:
%
%    fname      NetCDF file name (character string)
%    aname      Attribute name (character string)
%    vname      Variable name (character string; optional)
%
% On Output:
%
%    avalue     Attribute value (numeric or character string)
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%  Initialize.

avalue=[];

if (nargin < 3),
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

[ncid]=mexnc('open',fname,'nc_nowrite');
if (ncid < 0),
  disp('  ');
  error(['NC_GETATT: open - unable to open file: ', fname]);
  return
end

%---------------------------------------------------------------------------
%  Read requested variable attribute.
%---------------------------------------------------------------------------

if (got_var),

%  Get variable ID.

  [varid]=mexnc('inq_varid',ncid,vname);
  if (varid < 0),
    [status]=mexnc('close',ncid);
    disp('  ');
    error(['NC_GETATT: inq_varid - cannot find variable: ',vname]);
    return
  end,

%  Inquire number of variable attributes.

  [nvatts,status]=mexnc('inq_varnatts',ncid,varid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_GETATT: inq_varnatts - unable to inquire number of variable ' ...
	   'attributes: ',vname]);
  end,

%  Inquire if variable attribute exist.

  found=0;
  
  for i=0:nvatts-1
    [attnam,status]=mexnc('inq_attname',ncid,varid,i);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_GETATT: inq_attname: error while inquiring attribute: ' ...
	     num2str(i)]);
    end,
    lstr=length(attnam);
    if (strmatch(aname(1:latt),attnam(1:lstr),'exact')),
      found=1;
      break
    end,
  end,

%  Read in requested attribute.
  
  if (found),
    [atype,status]=mexnc('inq_atttype',ncid,varid,aname);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_GETATT: inq_atttype - unable to inquire datatype for ' ...
	     'attribute: ',aname]);
    end,

    switch (atype)
      case (ncchar)
        [avalue,status]=mexnc('get_att_text',ncid,varid,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_text - unable to read attribute ',...
                 '"',aname,'"in variable: ',vname,'.']);
        end,
      case (ncint)
        [avalue,status]=mexnc('get_att_int', ncid,varid,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_int - unable to read attribute ',...
                 '"',aname,'"in variable: ',vname,'.']);
        end,
      case (ncfloat)
        [avalue,status]=mexnc('get_att_float',ncid,varid,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_float - unable to read attribute ',...
                 '"',aname,'"in variable: ',vname,'.']);
        end,
      case (ncdouble)
        [avalue,status]=mexnc('get_att_double',ncid,varid,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_double - unable to read attribute ',...
                 '"',aname,'"in variable: ',vname,'.']);
        end,
    end,
  end,

%---------------------------------------------------------------------------
%  Read requested global attribute.
%---------------------------------------------------------------------------
   
else,

%  Inquire number of global attributes.

  [natts,status]=mexnc('inq_natts',ncid);
  if (status < 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_GETATT: inq_natts - unable to inquire number of global ' ...
	   'attributes: ',fname]);
  end,
  
%  Inquire if requested global attribute exist.

  found=0;
  
  for i=0:natts-1
    [attnam,status]=mexnc('inq_attname',ncid,ncglobal,i);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_GETATT: inq_attname: error while inquiring attribute: ' ...
	     num2str(i)]);
      return
    end,
    lstr=length(attnam);
    if (strmatch(aname(1:latt),attnam(1:lstr),'exact')),
      found=1;
      break
    end,
  end,
  
%  Read in requested attribute.

  if (found),
    [atype,status]=mexnc('inq_atttype',ncid,ncglobal,aname);
    if (status < 0),
      disp('  ');
      disp(mexnc('strerror',status));
      error(['NC_GETATT: inq_atttype - unable to inquire datatype for ' ...
	     'attribute: ',aname]);
    end,

    switch (atype)
      case (ncchar)
        [avalue,status]=mexnc('get_att_text',ncid,ncglobal,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_text - unable to read global ',...
                 'attribute: "',aname,'".']);
        end,
      case (ncint)
        [avalue,status]=mexnc('get_att_int', ncid,ncglobal,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_int - unable to read global ',...
                 'attribute: "',aname,'".']);
        end,
      case (ncfloat)
        [avalue,status]=mexnc('get_att_float',ncid,ncglobal,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_float - unable to read global ',...
                 'attribute: "',aname,'".']);
        end,
      case (ncdouble)
        [avalue,status]=mexnc('get_att_double',ncid,ncglobal,aname);
        if (status < 0),
          disp('  ');
          disp(mexnc('strerror',status));
          error(['NC_GETATT: get_att_double - unable to read global ',...
                 'attribute: "',aname,'".']);
        end,
    end,
  end,
  
end,

%  Close NetCDF file.

[cstatus]=mexnc('ncclose',ncid);
if (cstatus < 0),
  disp('  ');
  disp(mexnc('strerror',status));
  error(['NC_GETATT: ncclose - unable to close NetCDF file: ', fname]);
end

return
