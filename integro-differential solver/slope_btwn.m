%Author: Russell J. Phelan 
%Date: 1-28-16

function slope = slope_btwn(f,next,prev)
%returns slope btwn ith and jth index of function object. 

dx = f(2,next)-f(2,prev);

slope = (f(1,next)-f(1,prev))/dx;

end

