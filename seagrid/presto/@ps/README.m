README -- Documentation for the "ps" class.
 
The "ps" class in the "Presto" toolbox replaces the "px"
class of the "Proxy" toolbox.  The new class is much
simpler in design, easier to maintain, and extensible.

Two features make "ps" objects attractive:

	(1)	Their data can be extended arbitrarily, using simple
		Matlab "." syntax for assigning and referencing the
		fields.  Whenever such a field-name is thename of
		an associated graphics property, the property itself
		is targetted.  Unlike Matlab's "set/get" scheme, "ps"
		subscripting can be nested to any depth that makes
		sense, so that specific portions of a property or
		field can be accessed directly.  Thus, "myPS.color(2)"
		could refer to the second element of color in (say)
		a figure attached to myPS, a "ps" object.
		
	(2)	Classes can be derived from "ps" without restricting
		access to any of the user's data.  This behavior is
		unlike that seen in regular Matlab objects, whose
		data is cosidered "private".
		
The "ps" class is designed for handling events generated
by graphical user interfaces, such as dialogs.  Menus
and controls are easily created and connected to particular
handlers (methods).  Events are typically funneled into
a method called "doevent", whose argument is usually the
name of the callback, as in "psevent ResizeFcn".

The "ps_test.m" demonstration program shows hows to install
menus and controls in a figure.

Needed:

1.	A "bind()" routine, to be applied after the "ps" object
	is created, so that the correct object is stored in the
	"UserData" of the associated handle.
	
2.	Inheritance in "PST".

% svn $Id$
%=======================================================================
% Copyright (C) 1999 Dr. Charles R. Denham, ZYDECO.
%  All Rights Reserved.
%   Disclosure without explicit written consent from the
%    copyright owner does not constitute publication.
%=======================================================================
 
% Version of 09-Dec-1999 01:55:25.
% Updated    09-Dec-1999 01:55:25.
