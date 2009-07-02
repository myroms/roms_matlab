% Simple test of Mexrect
% 

% svn $Id$

x=-pi:.2:pi;
n=length(x);
y=tanh(x);
xbound=[ -pi x fliplr(x)];
ybound=[ -.05 y fliplr(y+1)];
zbound=zeros(size(xbound));
zbound([2 n+1 n+2 n*2+1])=1;
ind=find(zbound);
grid_outline=[xbound(:) ybound(:) zbound(:)];

imax = 1000;
xb = grid_outline(:,1); 
yb = grid_outline(:,2); 
icorner = grid_outline(:,3);

z = xb + i*yb;
n = zeros(4,1);
corner_find = find(icorner==1);
n = corner_find;


np = length(xb);


%
% map physical boundary to a rectangle
itmax = 80;
errmax = 1e-5;
error_within_bounds = 0;
zold = z;
for k = 1:itmax
	
	z = mexrect ( z, np, n(1), n(2), n(3), n(4) );

	%
	% calculate departure of contour from rectangle
	error = 0.0;
	error = error + sum ( abs ( real (z(1:n(1)) - z(1)) ) );
	error = error + sum ( abs ( imag (z(n(1)+1:n(2)) - z(n(1)+1)) ) );
	error = error + sum ( abs ( real (z(n(2)+1:n(3)) - z(n(2)+1)) ) );
	error = error + sum ( abs ( imag (z(n(3)+1:n(4)) - z(n(3)+1)) ) );

	error = error / (imag(z(n(4))) * 2 + 2 );

	reg_error(k) = error;
%	disp ( sprintf ( 'regularity error in mapped contour at iteration %4.0f is %f\n', k, error ) );

	if ( abs(error) < errmax )
		error_within_bounds = 1;
		break;
	end

end

figure;
semilogy([1:k],reg_error);
h = line('xdata',[1:k],'ydata',reg_error,'marker','*','color','r','linestyle','none');

title_str = sprintf ( 'Test of Mexrect:  iterative error should drop below %g.', errmax );
title(title_str);
xlabel ( 'Iteration' );
ylabel ( 'Error' );
set ( gca, 'xtick', [1:1:k] );

if ( error_within_bounds )
	;
else
	disp ( 'warning:  failed to converge in %4.0f iterations\n', itmax );
end

