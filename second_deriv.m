%Author: Russell J. Phelan 
%Date: 10-5-15

function yPrime = second_deriv(t)
 %a very sketchy second derivative function, used for putting an upper
 %bound on the error in simpson's rule. 
 
 h = .0001;
 
 yPrime = (deriv_at_point(t+h)-deriv_at_point(t))/h;
end

