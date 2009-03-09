
% SW_NEW    What's new in this version of seawater.
%
% 94/11/15 release 1.2d
% **********************
% sw_bfrq.m   Now also returns potential vorticity.
%             Thanks to Greg Johnson (gjohnson@pmel.noaa.gov)
%
% sw_gvel.m   OMEGA=7.29e-5 changed to OMEGA=7.292e-5 to be
%             consistent with sw_f.m
%
%             IMPORTANT CHANGE: The usage of the following 
%             routines has changed!
%
% sw_alpha.m |    All these routines expect (S,T,P) to
% sw_beta.m  |--  be passed instead of (S,PTMP,P) as in 
% sw_aonb.m  |    previous releases of seawater.
%                 Fast execution can still be obtained by passing
%                 ptmp with a string flag 'ptmp' see help.
%
% 94/10/19 release 1.2c
% **********************
% Added routine sw_new.m to inform of updates and new features.
% sw_bfrq.m   Fixed bug where LAT = [] was needed as argument when
%             no latitude values are being passed.
%             Now pass PRESSURE instead of DEPTH -> more consistent
%             though only a negligible change is answers.
%
% sw_info.m   Updated to include a registration section.
%             Noted that software is FREE.  
%             Noted best email address is seawater@ml.csiro.au
%             Requests for Report also via email to library@ml.csiro.au
%
% 94/10/12 release 1.2b
% ********************
% First official release and announcement on the networks.
%

more on
help sw_new
more off
