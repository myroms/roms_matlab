function editscope(grid_file, coast_file);

% EDITSCOPE: interactive adjoint sensitivity scope mask editing for ROMS
%
% editscope(grid_file, coast_file)
%
% GUI for manual editing of ROMS adjoint sensitivity scope mask on
% RHO-points.  To accelerate the processing, the Land/Sea mask is
% edited in (I,J) grid coordinates. If the (I,J) coordinates are not
% provided, it computes and writes them into file. If called without,
% one or both arguments, it will prompt for the needed file name(s).
% If the coastline data is in grid_file, it will read it and convert
% it to (I,J) fractional coordinates.
%
% On Input:
%
%    grid_file     ROMS Grid NetCDF file name containing the grid
%                  and mask arrays (string)
%
%    coast_file    Matlab file name containing the coastline
%                  (lon,lat) or (I,J) fractional coordinates
%                  (string; optional)
%
% Mouse shortcuts:
%
%    double click ==> Zoom in
%    right  click ==> Zoom out
%    middle click ==> change editing mode
%
%    Calls: READ_MASK, WRITE_MASK and UVP_MASKS functions.
%           BUTTON, RADIOBOX, TEXTBOX, AXISSCROLL,
%

% svn $Id$
%=========================================================================%
%  Copyright (c) 2002-2013 The ROMS/TOMS Group                            %
%    Licensed under a MIT/X style license                                 %
%    See License_ROMS.txt                           Hernan G. Arango      %
%=========================================================================%

% Define and initialize persistent variables.

persistent changed Rscope rlon rlat bath scope hplot hcst ha Lp Mp mx my
persistent mfile fig zooming xx yy xl yl xcst ycst

% Single global variable to pass info to/from the callback routines.

global GUI

% Set colormap: first two entries - land, second pair - sea

%CMAP=[0 1 0;.5 1 0;0 .7 1;.3 0 1];
CMAP=[.5 1 0;1 1 0;0 0 .7;0 0 1];

% Set coastline line color and width.

LineColor='k';
LineWidth=1;

% Check input arguments.

if (nargin < 1 | isempty(grid_file)),
  grid_file='*.nc';
end,
if (nargin < 2 | isempty(coast_file)),
  coast_file='*.mat';
end,
got_coast=0;

FIGURE_NAME='Scope Mask Editor';

%==========================================================================
% The EDITSCOPE function is also used as a CallBack for some of the
% uicontrols,  so we need to figure out, what's required by looking
% at the 'grid_file' parameter.
%==========================================================================

switch lower(grid_file),

%--------------------------------------------------------------------------
% Zoom-in.
%--------------------------------------------------------------------------

  case 'zoomin',

    disable_click;
    editscope move;
    zooming=1;
    waitforbuttonpress
    xx0=xx; yy0=yy;                         % save current pointer position
    rbbox;                                  % track rubberband rectangle
    xx1=xx; yy1=yy;                         % new pointer position
    if ((xx0 ~= xx1) & (yy0 ~= yy1)),       % trim limits and set zoom
      xx0(xx0<0)=0;  xx0(xx0>xl(2))=xl(2);
      xx1(xx1<0)=0;  xx1(xx1>xl(2))=xl(2);
      yy0(yy0<0)=0;  yy0(yy0>yl(2))=yl(2);
      yy1(yy1<0)=0;  yy1(yy1>yl(2))=yl(2);
      xlim([min([xx0 xx1]), max([xx0 xx1])]);
      ylim([min([yy0 yy1]), max([yy0 yy1])]);
      axisscroll;
    end,
    enable_click;
    zooming=0;

%--------------------------------------------------------------------------
% Zoom-out.
%--------------------------------------------------------------------------

  case 'zoomout',

    xlim(xl);
    ylim(yl);
    axisscroll;

%--------------------------------------------------------------------------
% Edit adjoint sensitivity scope mask.
%--------------------------------------------------------------------------

  case 'click',

    button=get(gcf, 'SelectionType');
    if (strcmp(button,'alt')),              % zoom out on right click
      editscope zoomout;
      return;
    end,
    if (strcmp(button,'open')),             % zoom in on double click
      editscope zoomin;
      return;
    end,
    if (strcmp(button,'extend')),           % cycle modes on middle click
      m=mod(GUI.mode,3)+1;
      eval(get(GUI.mode_h(m),'callback'));
      editscope zoomout;
      return;
    end,
    if (within(xx,xlim) & within(yy,ylim)), % left click within edit area
      disable_click;
      switch GUI.tool
        case 1,                             % point edit
          ix=floor(xx+1.5);
          iy=floor(yy+1.5);
          switch GUI.mode
            case 1,                         % toggle between off and on
              Rscope(ix,iy)=~Rscope(ix,iy);
            case 2,                         % off
              Rscope(ix,iy)=0;
            case 3,                         % on
              Rscope(ix,iy)=1;
	  end,
        case 2,                             % area edit
          xx0=xx; yy0=yy;                   % save current pointer position
          rbbox;                            % track rubberband rectangle
          xx1=xx; yy1=yy;                   % new pointer position
          idsel=find((mx-xx0-1).*(mx-xx1-1)<=0 & ...
                     (my-yy0-1).*(my-yy1-1)<=0);% indicies within rectangle
          switch GUI.mode
            case 1,
              Rscope(idsel)=~Rscope(idsel);      % toggle between off and on
            case 2,
              Rscope(idsel)=0;                   % off
            case 3,
              Rscope(idsel)=1;                   % on
          end,
      end,
      changed=1;
      enable_click;
      editscope refresh;
    end,

%---------------------------------------------------------------------------
% Update scope by changing colormap.
%---------------------------------------------------------------------------

  case 'refresh',

    cdata=Rscope*2+mod(mod(mx,2)+mod(my,2),2)+1;
    set(hplot,'cdata',cdata');
    nm=[FIGURE_NAME,' - ', mfile];
    if (changed),
      nm=[nm,'*'];
    end,
    set(fig,'name',nm);

%---------------------------------------------------------------------------
% Pointer movement: update pointer coordinates on menu.
%---------------------------------------------------------------------------

  case 'move',

    xy=get(gca,'currentpoint');
    xx=xy(1,1); yy=xy(1,2);
    if (within(xx,xlim) & within(yy,ylim)),
      s=sprintf('[%d,%d]',floor(xx+.5),floor(yy+.5));
      if (zooming),
        pointer('zoom');
      elseif (GUI.tool==1),
        pointer('point');
      elseif GUI.tool==2,
        pointer('rect');
      end,
    else,
      s='---';
      pointer('arrow');
    end,
    set(GUI.pos_h,'string',s);

%---------------------------------------------------------------------------
% Compute U and V scope masks.  Write out into GRID NetCDF file.
%---------------------------------------------------------------------------

  case 'save',

    [Uscope,Vscope]=uv_scope(Rscope);
    scope=Rscope;
    write_scope(mfile,Rscope,Uscope,Vscope);
    changed=0;
    editscope refresh;

%---------------------------------------------------------------------------
% Undo changes: restore last saved scope mask.
%---------------------------------------------------------------------------

  case 'undo',

    Rscope=scope;
    changed=0;
    editscope refresh;

%---------------------------------------------------------------------------
% Done: inquire to save scope mask changes.
%---------------------------------------------------------------------------

  case 'exit',

    if (~changed),
      delete(gcf);
    else
      res=questdlg('The scope mask has been changed. Save?',FIGURE_NAME);
      switch res
        case 'Yes',
          editscope save;
          disp('Scope mask has been saved');
          delete(gcf);
        case 'No',
          disp('Scope mask has NOT been saved');
          delete(gcf);
      end,
    end,

%---------------------------------------------------------------------------
% Help.
%---------------------------------------------------------------------------

  case 'help',

    show_help;

%---------------------------------------------------------------------------
% Initialize: read scope mask and coastline data.
%---------------------------------------------------------------------------

  otherwise,

% Kill all windows.

    delete(findobj(0,'tag','ScopeEditor'));

% If appropriate, inquire input file names.

    if (any(grid_file=='*')),
      [fn,pth]=uigetfile(grid_file,'Select ROMS grid file...');
      if (~fn),
        return,
      end;
      grid_file=[pth,fn];
    end,

% Read in grid data.

   mfile=grid_file;
   [spherical,rlon,rlat,bath,scope]=read_scope(grid_file);
   Rscope=scope;
   [Lp,Mp]=size(scope);
   [mx,my]=ndgrid(1:Lp,1:Mp);

% Check if grid NetCDF file has coastline data.

   S=nc_vnames(grid_file);
   vnames={S.Variables.Name};
   got_Clon=any(strcmp(vnames,'lon_coast'));
   got_Clat=any(strcmp(vnames,'lat_coast'));

   if (~(got_Clon & got_Clat) & spherical),
     if (any(coast_file=='*')),
       [fn,pth]=uigetfile(coast_file,'Select Coastline file...');
       if (~fn),
%        return,
         got_coast=0;
       else,
         got_coast=1;
       end;
       coast_file=[pth,fn];
     end,
   end,

% Read in coastline data. If appropriate, compute coastline (I,J)
% grid indices.

   if (got_Clon & got_Clat),
     [C]=ijcoast(grid_file);
     xcst=C.Icst;
     ycst=C.Jcst;
     clear C lat lon;
   else,
     if (got_coast),
       load(coast_file);
       if (exist('C')),
         xcst=C.Icst;
         ycst=C.Jcst;
         clear C;
       elseif (exist('lon') & exist('lat')),
         [C]=ijcoast(grid_file,coast_file);
         xcst=C.Icst;
         ycst=C.Jcst;
         clear C lat lon;
       else,
         error('Coast file should contain "lon" and "lat" vectors');
       end,
     end,
   end,

% Initialize the window.

   fig=figure('NumberTitle','off',...
              'tag','ScopeEditor',...
              'DoubleBuffer','on',...
              'backingstore','off',...
              'menubar','none');
   cdata=Rscope*2+mod(mod(mx,2)+mod(my,2),2)+1;
   gx=mx(:,1)-1;
   gy=my(1,:)-1;
   axes('position',[0.0554 0.1024 0.7518 0.8405]);
   hplot=image(gx,gy,cdata','cdatamapping','direct','erasemode','normal');
   set(gca,'YDir','normal',...
           'drawmode','fast',...
           'layer','top',...
           'tickdir','out');
   colormap(CMAP);
   hold on;
   hline=plot(xcst,ycst,LineColor);
   set(hline,'LineWidth',LineWidth);
   xl=xlim; yl=ylim;
   changed=0;
   setgui;
   editscope refresh;

end,

return


function setgui

%---------------------------------------------------------------------------
% set-up scope mask editing menu bottons.
%---------------------------------------------------------------------------

textbox([.85 .85 .145 .14], ....
        '(i,j)',{'0,0'},'pos');

radiobox([.85 .65 .145 .19], ...
         'Edit Mode',{'Toggle Scope','off','on'},'mode');

radiobox([.85 .5 .145 .14], ...
         'Edit Tool',{'Point edit','Area edit'},'tool');

button([.85 .4 .145 .05], ...
       'Zoom In',[mfilename ' zoomin']);

button([.85 .35 .145 .05], ...
       'Zoom Out',[mfilename ' zoomout']);

button([.85 .25 .145 .05], ...
       'Undo',[mfilename ' undo']);

button([.85 .2 .145 .05], ...
       'Save',[mfilename ' save']);

button([.85 .1 .145 .05], ...
       'Exit',[mfilename ' exit']);

axisscroll('r')
axisscroll('t')
set(gcf,'WindowButtonMotionFcn',[mfilename ' move;'],...
        'CloseRequestFcn',[mfilename ' exit;'],...
        'interruptible','on');
enable_click;

return

function disable_click

%---------------------------------------------------------------------------
% Disable pointer clicking on current figure.
%---------------------------------------------------------------------------

set(gcf,'WindowButtonDownFcn','');

return

function enable_click

%---------------------------------------------------------------------------
% Enable pointer clicking on current figure.
%---------------------------------------------------------------------------

set(gcf,'WindowButtonDownFcn',[mfilename ' click;']);

return

function r=within(a,b)


%---------------------------------------------------------------------------
% Check if 'a' is within the range of 'b'.
%---------------------------------------------------------------------------

r=(a>=b(1) & a<= b(2));

return
