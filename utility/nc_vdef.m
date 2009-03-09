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
%                  Var.name  => variable name.
%                  Var.type  => variable type.
%                  Var.dimid => variable dimension IDs.
%                  Var.long  => variable "long-name" attribute.
%                  Var.opt_T => variable "option_T" attribute.
%                  Var.opt_F => variable "option_F" attribute.
%                  Var.opt_0 => variable "option_0" attribute.
%                  Var.opt_1 => variable "option_1" attribute.
%                  Var.opt_2 => variable "option_2" attribute.
%                  Var.opt_3 => variable "option_3" attribute.
%                  Var.opt_4 => variable "option_4" attribute.
%                  Var.opt_5 => variable "option_5" attribute.
%                  Var.opt_6 => variable "option_6" attribute.
%                  Var.opt_7 => variable "option_7" attribute.
%                  Var.units => variable "units" attribute.
%                  Var.fill  => variable "_FillValue" attribute.
%                  Var.miss  => variable "missing_value" attribute.
%                  Var.offset=> variable "add_offset" attribute.
%                  Var.min   => variable "valid_min" attribute.
%                  Var.max   => variable "valid_max" attribute.
%                  Var.minus => variable "negative" attribute.
%                  Var.plus  => variable "positive" attribute.
%                  Var.field => variable "field" attribute.
%                  Var.pos   => variable "positions" attribute.
%                  Var.time  => variable "time" attribute.
%                  Var.coord => variable "coordinates" attribute.
%                  Var.urot  => variable "u-rotation" attribute.
%                  Var.vrot  => variable "v-rotation" attribute.
%                  Var.left  => variable "left" attribute.
%                  Var.right => variable "right" attribute.
%                  Var.top   => variable "top" attribute.
%                  Var.bottom=> variable "bottom" attribute.
%                  Var.up    => variable "up" attribute.
%                  Var.down  => variable "down" attribute.
%
% On Output:
%
%    varid       Variable ID.
%    status      Error flag.
%

% svn $Id$
%===========================================================================%
%  Copyright (c) 2002-2009 The ROMS/TOMS Group                              %
%    Licensed under a MIT/X style license                                   %
%    See License_ROMS.txt                           Hernan G. Arango        %
%===========================================================================%


%---------------------------------------------------------------------------
%  Get some NetCDF parameters.
%---------------------------------------------------------------------------

[ncdouble]=mexnc('parameter','nc_double');
[ncfloat ]=mexnc('parameter','nc_float');
[ncchar  ]=mexnc('parameter','nc_char');

% Set variable type (default: floating point, single precision)

if (isfield(Var,'type')),
  vartyp=Var.type;
else,
  vartyp=ncfloat;
end,

%---------------------------------------------------------------------------
%  Define requested variable.
%---------------------------------------------------------------------------

% Define variable.

if (isfield(Var,'name') & isfield(Var,'dimid')),
  vdid=Var.dimid;
  nvdim=length(vdid);
  [varid]=mexnc('ncvardef',ncid,Var.name,vartyp,nvdim,vdid);
  if (varid == -1),
    error(['NC_VDEF: ncvardef - unable to define variable: ',Var.name]);
    return,
  end,
end,

% Set "long-name" attribute.

if (isfield(Var,'long')),
  text=Var.long;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'long_name',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Var.name,':long_name.']);
      return,
    end,
  end,
end,

% Set "option_T" attribute.

if (isfield(Var,'opt_T')),
  text=Var.opt_T;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_T',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_T.']);
      return,
    end,
  end,
end,

% Set "option_F" attribute.

if (isfield(Var,'opt_F')),
  text=Var.opt_F;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_F',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_T.']);
      return,
    end,
  end,
end,

% Set "option_0" attribute.

if (isfield(Var,'opt_0')),
  text=Var.opt_0;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_0',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_0.']);
      return,
    end,
  end,
end,

% Set "option_1" attribute.

if (isfield(Var,'opt_1')),
  text=Var.opt_1;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_1',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_1.']);
      return,
    end,
  end,
end,

% Set "option_2" attribute.

if (isfield(Var,'opt_2')),
  text=Var.opt_2;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_2',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_2.']);
      return,
    end,
  end,
end,

% Set "option_3" attribute.

if (isfield(Var,'opt_3')),
  text=Var.opt_3;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_3',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_3.']);
      return,
    end,
  end,
end,

% Set "option_4" attribute.

if (isfield(Var,'opt_4')),
  text=Var.opt_4;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_4',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_4.']);
      return,
    end,
  end,
end,

% Set "option_5" attribute.

if (isfield(Var,'opt_5')),
  text=Var.opt_5;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_5',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_5.']);
      return,
    end,
  end,
end,

% Set "option_6" attribute.

if (isfield(Var,'opt_6')),
  text=Var.opt_6;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_6',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_6.']);
      return,
    end,
  end,
end,

% Set "option_7" attribute.

if (isfield(Var,'opt_7')),
  text=Var.opt_7;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'option_7',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_7.']);
      return,
    end,
  end,
end,

% Set "units" attribute.

if (isfield(Var,'units')),
  text=Var.units;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'units',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Var.name,':units.']);
      return,
    end,
  end,
end,

% Set "add_offset" attribute.

if (isfield(Var,'offset')),
  [status]=mexnc('ncattput',ncid,varid,'add_offset',vartyp,1,Var.offset);
  if (status == -1),
    error(['NC_VDEF: ncattput - unable to define attribute: ',...
           Vname.name,':valid_offset']);
    return,
  end,
end,

% Set "valid_min" attribute.

if (isfield(Var,'min')),
  [status]=mexnc('ncattput',ncid,varid,'valid_min',vartyp,1,Var.min);
  if (status == -1),
    error(['NC_VDEF: ncattput - unable to define attribute: ',...
           Vname.name,':valid_min']);
    return,
  end,
end,

% Set "valid_max" attribute.

if (isfield(Var,'max')),
  [status]=mexnc('ncattput',ncid,varid,'valid_max',vartyp,1,Var.max);
  if (status == -1),
    error(['NC_VDEF: ncattput - unable to define attribute: ',...
           Vname.name,':valid_max.']);
    return,
  end,
end,

% Set "positive" attribute.

if (isfield(Var,'plus')),
  text=Var.plus;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'positive',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_T.']);
      return,
    end,
  end,
end,

% Set "negative" attribute.

if (isfield(Var,'minus')),
  text=Var.minus;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'negative',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':option_T.']);
      return,
    end,
  end,
end,

% Set "_FillValue" attribute.

if (isfield(Var,'fill')),
  [status]=mexnc('ncattput',ncid,varid,'_FillValue',vartyp,1,Var.fill);
  if (status == -1),
    error(['NC_VDEF: ncattput - unable to define attribute: ',...
           Vname.name,':_FillValue']);
    return,
  end,
end,

% Set "missing_value" attribute.

if (isfield(Var,'miss')),
  [status]=mexnc('ncattput',ncid,varid,'missing_value',vartyp,1,Var.miss);
  if (status == -1),
    error(['NC_VDEF: ncattput - unable to define attribute: ',...
           Vname.name,':missing_value']);
    return,
  end,
end,

% Set "positions" attribute.

if (isfield(Var,'pos')),
  text=Var.pos;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'positions',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':positions']);
    end,
  end,
end,

% Set "time" attribute.

if (isfield(Var,'time')),
  text=Var.time;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'time',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':time']);
    end,
  end,
end,

% Set "coordinates" attribute.

if (isfield(Var,'coord')),
  text=Var.coord;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'coordinates',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':coordinates']);
    end,
  end,
end,

% Set "u-rotation" attribute.

if (isfield(Var,'urot')),
  text=Var.pos;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'rotation1',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':rotation1']);
    end,
  end,
end,

% Set "v-rotation" attribute.

if (isfield(Var,'vrot')),
  text=Var.pos;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'rotation2',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':rotation2']);
    end,
  end,
end,

% Set "left" attribute.

if (isfield(Var,'left')),
  text=Var.left;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'left',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':left.']);
      return,
    end,
  end,
end,

% Set "right" attribute.

if (isfield(Var,'right')),
  text=Var.right;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'right',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':right.']);
      return,
    end,
  end,
end,

% Set "top" attribute.

if (isfield(Var,'top')),
  text=Var.top;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'top',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':top.']);
      return,
    end,
  end,
end,

% Set "bottom" attribute.

if (isfield(Var,'bottom')),
  text=Var.bottom;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'bottom',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':bottom.']);
      return,
    end,
  end,
end,


% Set "up" attribute.

if (isfield(Var,'up')),
  text=Var.up;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'up',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':up.']);
      return,
    end,
  end,
end,

% Set "down" attribute.

if (isfield(Var,'down')),
  text=Var.down;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'down',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':down.']);
      return,
    end,
  end,
end,

% Set "field" attribute.

if (isfield(Var,'field')),
  text=Var.field;
  lstr=length(text);
  if (lstr > 0),
    [status]=mexnc('ncattput',ncid,varid,'field',ncchar,lstr,text);
    if (status == -1),
      error(['NC_VDEF: ncattput - unable to define attribute: ',...
             Vname.name,':field.']);
      return,
    end,
  end,
end,

return
