%Author: Russell J. Phelan 
%Date: 1-28-16

function slope = slope_btwn(f,curr_t,prev_t)
%returns slope btwn ith and jth index of function object. 

dx = f(2,curr_t)-f(2,prev_t);

slope = (f(1,curr_t)-f(1,prev_t))/dx;

end

