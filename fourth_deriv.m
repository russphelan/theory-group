%Author: Russell J. Phelan 
%Date: 10-5-15

function yPrime = fourth_deriv(t)
%a very sketchy fourth derivative function, used for putting an upper
%bound on the error in simpson's rule. 

h = .01;

yPrime = (third_deriv(t+h)-third_deriv(t))/h;
end