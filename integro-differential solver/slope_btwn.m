%Author: Russell J. Phelan 
%Date: 8/12/16

function slope = slope_btwn(f,next,prev)
%returns slope at point between next and prev

dx = f(2,prev+1)-f(2,prev);

slope = (f(1,next)-f(1,prev))/2/dx;

end

