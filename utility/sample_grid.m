function [Istr,Iend,Jstr,Jend] = sample_grid(XP,YP,XT,YT,varargin)

%
% SAMPLE_GRID:  Gets parent grid tight indices containing target grid
%
% [Istr,Iend,Jstr,Jend] = sample_grid(XP,YP,XT,YT,offset,plt);
%
% This function computes the parent grid indices range of the polygon
% that tightly contains the target grid.  This is done to sample the
% parent grid to accelerate the interpolation of fields for the target
% grid.
%
% On Input:
%
%    XP            Parent grid X-coordinates (2D array)
%    YP            Parent grid Y-coordinates (2D array)
%    XT            Target grid X-coordinates (2D array)
%    YT            Target grid Y-coordinates (2D array)
%    offset        Number of extra points to used to sample the
%                    parent grid so is large enough to contain
%                    the target grid  (default 5)
%    plt           Switch to plot parent and target grids
%                    (default false)
%
% On Output:
%
%    Istr          Parent grid starting I-index for sampling
%    Iend          Parent grid ending   I-index for sampling
%    Jstr          Parent grid starting J-index for sampling
%    Jend          Parent grid ending   J-index for sampling
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2012 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

%  Set optional arguments.

offset = 5;
plt = false;

switch numel(varargin)
  case 1
    offset = varargin{1};
  case 2
    offset = varargin{1};
    plt = varargin{2};
end

%--------------------------------------------------------------------------
%  Extract parent grid polygon.
%--------------------------------------------------------------------------

[ImP,JmP] = size(XP);
IstrP = 1;
IendP = ImP;
JstrP = 1;
JendP = JmP;

XboxP = [squeeze(XP(IstrP:IendP,JstrP));                              ...
         squeeze(XP(IendP,JstrP+1:JendP))';                           ...
         squeeze(flipud(XP(IstrP:IendP-1,JendP)));                    ...
         squeeze(fliplr(XP(IstrP,JstrP:JendP-1)))'];

YboxP = [squeeze(YP(IstrP:IendP,JstrP));                              ...
         squeeze(YP(IendP,JstrP+1:JendP))';                           ...
         squeeze(flipud(YP(IstrP:IendP-1,JendP)));                    ...
         squeeze(fliplr(YP(IstrP,JstrP:JendP-1)))'];

%--------------------------------------------------------------------------
%  Extract target grid polygon.
%--------------------------------------------------------------------------

[ImT,JmT] = size(XT);
IstrT = 1;
IendT = ImT;
JstrT = 1;
JendT = JmT;

XboxT = [squeeze(XT(IstrT:IendT,JstrT));                              ...
         squeeze(XT(IendT,JstrT+1:JendT))';                           ...
         squeeze(flipud(XT(IstrT:IendT-1,JendT)));                    ...
         squeeze(fliplr(XT(IstrT,JstrT:JendT-1)))'];

YboxT = [squeeze(YT(IstrT:IendT,JstrT));                              ...
         squeeze(YT(IendT,JstrT+1:JendT))';                           ...
         squeeze(flipud(YT(IstrT:IendT-1,JendT)));                    ...
         squeeze(fliplr(YT(IstrT,JstrT:JendT-1)))'];

%--------------------------------------------------------------------------
%  Find parent indices containing the target grid.
%--------------------------------------------------------------------------

[in,on] = inpolygon(XP,YP,XboxT,YboxT);

%--------------------------------------------------------------------------
%  Set parent grid sampling indices.
%--------------------------------------------------------------------------

[J,I] = meshgrid(1:1:JmP, 1:1:ImP);

I(~in) = NaN;
J(~in) = NaN;

Istr = min(I(:))-offset;
if (isnan(Istr) || Istr < 1),
  Istr = 1;
end

Iend = max(I(:))+offset;
if (isnan(Iend) || Iend > ImP),
  Iend = ImP;
end

Jstr = min(J(:))-offset;
if (isnan(Jstr) || Jstr < 1),
  Jstr = JmP;
end

Jend = max(J(:))+offset;
if (isnan(Jend) || Jend > JmP),
  Jend = JmP;
end

%--------------------------------------------------------------------------
%  If requested, plot parent and target grids.
%--------------------------------------------------------------------------

if (plt),

  X = XP(Istr:1:Iend,Jstr:1:Jend);
  Y = YP(Istr:1:Iend,Jstr:1:Jend);

  figure;
  plot(XboxP,YboxP,'r.',XboxT,YboxT,'b.');
  title(['Original Parent and Taget Grids']);
  
  figure;
  pcolor(X,Y,ones(size(X)));
  hold on;
  plot(XboxT,YboxT,'b.');
  title(['Sampled Parent and Taget Grids']);
  hold off

end,

return
