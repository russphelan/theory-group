%Author: Russell J. Phelan 
%Date: 10-5-15

function yPrime = deriv_at_point(t)
%the function referenced is f.m in this package. 

%takes a point, returns a rough derivative at that point

%no error is included; this is used to put an upper bound on the
%error in simpson's rule

h = .0001;

%derivative forward:
yPrime = (f(t+h)-f(t))/h;
end