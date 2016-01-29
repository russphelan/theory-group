%Author: Russell J. Phelan 
%Date: 1-28-16

function slope = slope_btwn(f,i,j)
%returns slope btwn ith and jth index of function object. 

dx = f(2,j)-f(2,i);

slope = (f(1,j)-f(1,i))/dx;

end

