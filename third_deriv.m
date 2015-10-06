%Author: Russell J. Phelan 
%Date: 10-5-15

function yPrime = third_deriv(t)
%a very sketchy third derivative function, used for putting an upper
%bound on the error in simpson's rule. 
 
 h = .001;
 
 yPrime = (second_deriv(t+h)-second_deriv(t))/h;
end

