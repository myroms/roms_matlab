function theResult = dohelp(self)

% dohelp -- Show SeaGrid help dialog.
%  dohelp(self) displays the "help" dialog
%   for self, a "seagrid" object.
 
% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 20-Dec-1999 14:29:59.
% Updated    10-Jan-2000 18:04:11.

theClass = class(self);
theMethods = sort(methods(theClass));

for i = length(theMethods):-1:1
	theName = theMethods{i};
	if ~any(findstr(theMethods{i}, 'help_'))
		theMethods(i) = [];
	end
end

theHelp = theMethods;

if (0)
	
	s = [];
	t = [];
	i = 0;
	while i < length(theHelp)/2-1
		i = i+1;
		s = setfield(s, theHelp{i}, help(theHelp{i}));
	end
	while i < length(theHelp)
		i = i+1;
		t = setfield(t, theHelp{i}, help(theHelp{i}));
	end
	s.More_Topics = t;

else
	
	s = [];
	i = 0;
	k = 0;
	while i < length(theHelp)
		k = k+1;
		s{k} = [];
		j = 0;
		while j < 10 & i < length(theHelp)
			i = i+1;
			j = j+1;
			s{k} = setfield(s{k}, theHelp{i}, help(theHelp{i}));
		end
	end
	for k = length(s)-1:-1:1
		s{k} = setfield(s{k}, ['More_Topics_' int2str(k)], s{k+1});
	end
	s = s{1};

end

result = guido(s, 'SeaGrid Help');

if nargout > 0
	theResult = result;
end
