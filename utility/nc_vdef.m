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
%                  Var.name     => name (string)
%                  Var.type     => type (number)
%                  Var.dimid    => dimension IDs (number)
%                  Var.long     => "long_name" attribute (string)
%                  Var.flag_str => "flag_values" attribute (string)
%                  Var.flag_num => "flag_values" attribute (number)
%                  Var.meaning  => "flag_meanings" attribute (string)
%                  Var.units    => "units" attribute (string)
%                  Var.calendar => "calendar" attribute (string) 
%                  Var.offset   => "add_offset" attribute (number)
%                  Var.cycle    => "cycle_length" attribute (number)
%                  Var.min      => "valid_min" attribute (number)
%                  Var.max      => "valid_max" attribute (number)
%                  Var.positive => "positive" attribute (string)
%                  Var.plus     => "positive_value" attribute (string)
%                  Var.minus    => "negative_value" attribute (string)
%                  Var.fill     => "_FillValue" attribute (number)
%                  Var.miss     => "missing_value" attribute (string)
%                  Var.stdname  => "standard_name" attribute (string)
%                  Var.formula  => "formula_terms" attribute (string
%                  Var.time     => "time" attribute (string)
%                  Var.pos      => "positions" attribute (string)
%                  Var.coord    => "coordinates" attribute (string)
%                  Var.urot     => "u-rotation" attribute (string)
%                  Var.vrot     => "v-rotation" attribute (string)
%                  Var.left     => "left" attribute (string)
%                  Var.right    => "right" attribute (string)
%                  Var.top      => "top" attribute (string)
%                  Var.bottom   => "bottom" attribute (string)
%                  Var.up       => "up" attribute (string)
%                  Var.down     => "down" attribute (string)
%                  Var.field    => "field" attribute (string)
%
% On Output:
%
%    varid       Variable ID.
%    status      Error flag.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2010 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%---------------------------------------------------------------------------
%  Get some NetCDF parameters.
%---------------------------------------------------------------------------

[ncdouble]=mexnc('parameter','nc_double');
[ncfloat ]=mexnc('parameter','nc_float');
[ncchar  ]=mexnc('parameter','nc_char');
[ncint   ]=mexnc('parameter','nc_int');

% Set variable type (default: floating point, single precision)

if (isfield(Var,'type')),
  vartyp=Var.type;
else,
  vartyp=ncfloat;
end,

%---------------------------------------------------------------------------
%  Define requested variable.
%---------------------------------------------------------------------------

if (isfield(Var,'name') & isfield(Var,'dimid')),
  vdid=Var.dimid;
  nvdim=length(vdid);
  [varid,status]=mexnc('def_var',ncid,Var.name,vartyp,nvdim,vdid);
  if (varid == -1 | status ~= 0),
    error(['NC_VDEF: DEF_VAR - unable to define variable: ',Var.name]);
    return,
  end,
end,

%---------------------------------------------------------------------------
% Write variable attributes.
%---------------------------------------------------------------------------

names=fieldnames(Var);
nfields=length(names);

for n=1:nfields,
  put_string=0;
  put_number=0;
  Aname=char(names(n));
  switch Aname
    case ('long')
      Vatt='long_name';
      text=getfield(Var,Aname);
      put_string=1;
    case ('flag_str')
      Vatt='flag_values';
      text=getfield(Var,Aname);
      put_string=1;
    case ('flag_num')
      Vatt='flag_values';
      value=getfield(Var,Aname)
      put_number=1;
    case ('meaning')
      Vatt='flag_meanings';
      text=getfield(Var,Aname);
      put_string=1;
    case ('units')
      Vatt='units';
      text=getfield(Var,Aname);
      put_string=1;
    case ('calendar')
      Vatt='calendar';
      text=getfield(Var,Aname);
      put_string=1;
    case ('offset')
      Vatt='add_offset';
      value=getfield(Var,Aname);
      put_number=1;
    case ('cycle')
      Vatt='cycle_length';
      value=getfield(Var,Aname);
      put_number=1;
    case ('min')
      Vatt='valid_min';
      value=getfield(Var,Aname);
      put_number=1;
    case ('max')
      Vatt='valid_max';
      value=getfield(Var,Aname);
      put_number=1;
    case ('positive')
      Vatt='positive';
      text=getfield(Var,Aname);
      put_string=1;
    case ('plus')
      Vatt='positive_value';
      text=getfield(Var,Aname);
      put_string=1;
    case ('minus')
      Vatt='negative_value';
      text=getfield(Var,Aname);
      put_string=1;
    case ('fill')
      Vatt='_FillValue';
      value=getfield(Var,Aname);
      put_number=1;
    case ('miss')
      Vatt='missing_value';
      value=getfield(Var,Aname);
      put_number=1;
    case ('stdname')
      Vatt='standard_name';
      text=getfield(Var,Aname);
      put_string=1;
    case ('formula')
      Vatt='formula_terms';
      text=getfield(Var,Aname);
      put_string=1;
    case ('time')
      Vatt='time';
      text=getfield(Var,Aname);
      put_string=1;
    case ('pos')
      Vatt='positions';
      text=getfield(Var,Aname);
      put_string=1;
    case ('coord')
      Vatt='coordinates';
      text=getfield(Var,Aname);
      put_string=1;
    case ('urot')
      Vatt='rotation1';
      text=getfield(Var,Aname);
      put_string=1;
    case ('vrot')
      Vatt='rotation2';
      text=getfield(Var,Aname);
      put_string=1;
    case ('left')
      Vatt='left';
      text=getfield(Var,Aname);
      put_string=1;
    case ('right')
      Vatt='right';
      text=getfield(Var,Aname);
      put_string=1;
    case ('top')
      Vatt='top';
      text=getfield(Var,Aname);
      put_string=1;
    case ('bottom')
      Vatt='bottom';
      text=getfield(Var,Aname);
      put_string=1;
    case ('up')
      Vatt='up';
      text=getfield(Var,Aname);
      put_string=1;
    case ('down')
      Vatt='down';
      text=getfield(Var,Aname);
      put_string=1;
    case ('field')
      Vatt='field';
      text=getfield(Var,Aname);
      put_string=1;
  end,
  if (put_string),
    lstr=length(text);
    [status]=mexnc('put_att_text',ncid,varid,Vatt,ncchar,lstr,text);
    if (status ~= 0),
      error(['NC_VDEF: PUT_ATT_TEXT - unable to define attribute: ',...
             Var.name,':',Vatt,'.']);
      return,
    end,
  end,
  if (put_number),
    nval=length(value);
    switch (vartyp)
      case (ncint)   
        value=int32(value);
        [status]=mexnc('put_att_int',   ncid,varid,Vatt,vartyp,nval,value);
        if (status == -1),
          error(['NC_VDEF: PUT_ATT_INT - unable to define attribute: ',...
                 Vname.name,':',Vatt,'.']);
          return,
        end,
      case (ncfloat)
        value=single(value);
	[status]=mexnc('put_att_float', ncid,varid,Vatt,vartyp,nval,value);
        if (status == -1),
          error(['NC_VDEF: PUT_ATT_FLOAT - unable to define attribute: ',...
                 Vname.name,':',Vatt,'.']);
          return,
        end,
      case (ncdouble)
        value=double(value);
	[status]=mexnc('put_att_double',ncid,varid,Vatt,vartyp,nval,value);
        if (status == -1),
          error(['NC_VDEF: PUT_ATT_DOUBLE - unable to define attribute: ',...
                 Vname.name,':',Vatt,'.']);
          return,
        end,
    end,
  end,
end,

return
