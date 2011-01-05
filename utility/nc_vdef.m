function [varid,status]=nc_vdef(ncid,Var);

%
% NC_VDEF:   Create a ROMS variable in a NetCDF file
%
% [varid,status]=nc_vdef(ncid,Var)
%
% This function defines a variable into NetCDF file.
%
% On Input:
%
%    ncid        NetCDF file ID (integer).
%    Var         Variable information (structure array):
%
%                  Var.name  => variable name (string)
%
%                  Var.type  => external data type (string or numeric):
%                               'byte'   or  nc_byte
%                               'char'   or  nc_char
%                               'short'  or  nc_short
%                               'int'    or  nc_int
%                               'float'  or  nc_float
%                               'double' or  nc_double
%                                     
%                  Var.dimid => dimension IDs (numeric), if not present
%                               or empthy [] the variable is a scalar
%                  
%                  Any other field in the structure, if any,
%                  is proccessed as a variable attribute. For
%                  example:
%
%                  Var.long_name  => "long_name" attribute (string)
%                  Var.add_offset => "add_offset" attribute (number)
%
%                  the values of the attribute can be numeric (scalar
%                  or vector) or characters (array or cell array).
%
% On Output:
%
%    varid       Variable ID.
%    status      Error flag.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2011 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%

% Check input structure.

if (~isfield(Var,'name')),
  disp(' ');
  error([ 'NC_VDEF - Cannot field ''name'' in structure array: Var']);
  return
end,

if (~isfield(Var,'type')),
  disp(' ');
  error([ 'NC_VDEF - Cannot field ''type'' in structure array: Var']);
  return
end,

if (~isfield(Var,'dimid')),
  disp(' ');
  error([ 'NC_VDEF - Cannot field ''dimid'' in structure array: Var']);
  return
end,

%---------------------------------------------------------------------------
%  Get some NetCDF parameters.
%---------------------------------------------------------------------------

[ncbyte  ]=mexnc('parameter','nc_byte');
[ncchar  ]=mexnc('parameter','nc_char');
[ncshort ]=mexnc('parameter','nc_short');
[ncint   ]=mexnc('parameter','nc_int');
[ncfloat ]=mexnc('parameter','nc_float');
[ncdouble]=mexnc('parameter','nc_double');

% Set variable external data type representation.

if (isfield(Var,'type')),
  type=getfield(Var,'type');
  if (ischar(type)),
    switch type
      case 'byte'
        vartyp=ncbyte;
      case 'char'
        vartyp=ncchar;
      case 'short'
        vartyp=ncshort;
      case 'int'
        vartyp=ncint;
      case 'float'
        vartyp=ncfloat;
      case 'double'
        vartyp=ncdouble;
    end,
  else
    vartyp=type;
  end,
else,
  error(['NC_VDEF: external data type field ''type'' is missing in ' ...
	 'structure: Var']);
  return,
end,


%---------------------------------------------------------------------------
%  Define requested variable.
%---------------------------------------------------------------------------

if (isfield(Var,'name') & isfield(Var,'dimid')),
  vdid=Var.dimid;
  nvdim=length(vdid);
  [varid,status]=mexnc('def_var',ncid,Var.name,vartyp,nvdim,vdid);
  if (varid == -1 | status ~= 0),
    disp('  ');
    disp(mexnc('strerror',status));
    error(['NC_VDEF: DEF_VAR - unable to define variable: ',Var.name]);
    return,
  end,
end,

%---------------------------------------------------------------------------
%  Add variable attributes.
%---------------------------------------------------------------------------

names=fieldnames(Var);
nfields=length(names);

for n=1:nfields,
  Aname=char(names(n));
  switch Aname
    case {'name', 'type', 'dimid'}
      put_attribute=0;
    otherwise,
      put_attribute=1;
      value=getfield(Var,Aname);    
  end,

% Define variable attributes.  

  if (put_attribute),

% Attribute value is character array or cell array.

    if (iscellstr(value) | ischar(value)),
      if (iscellstr(value)),
        value=char(value);            % need figure out this one
      end,
      lstr=length(value);
      [status]=mexnc('put_att_text',ncid,varid,Aname,ncchar,lstr,value);
      if (status ~= 0),
        disp('  ');
        disp(mexnc('strerror',status));
        error(['NC_VDEF: PUT_ATT_TEXT - unable to define attribute: ',...
               Var.name,':',Aname,'.']);
        return,
      end,

% Attribute value is numeric (scalar or vector). Check external data
% representation.

    else,

      nval=length(value);
      switch (vartyp)
        case (ncint)   
          value=int32(value);
          [status]=mexnc('put_att_int',   ncid,varid,Aname,vartyp,nval,value);
          if (status ~= 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_VDEF: PUT_ATT_INT - unable to define attribute: ',...
                   Vname.name,':',Aname,'.']);
            return,
          end,
        case (ncfloat)
          value=single(value);
	  [status]=mexnc('put_att_float', ncid,varid,Aname,vartyp,nval,value);
          if (status ~= 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_VDEF: PUT_ATT_FLOAT - unable to define attribute: ',...
                   Vname.name,':',Aname,'.']);
            return,
          end,
        case (ncdouble)
          value=double(value);
	  [status]=mexnc('put_att_double',ncid,varid,Aname,vartyp,nval,value);
          if (status ~= 0),
            disp('  ');
            disp(mexnc('strerror',status));
            error(['NC_VDEF: PUT_ATT_DOUBLE - unable to define attribute: ',...
                   Vname.name,':',Aname,'.']);
            return,
          end,
      end,
    end,
  end,
end,

return
