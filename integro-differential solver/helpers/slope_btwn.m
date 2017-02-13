%Author: Russell J. Phelan 
%Date: 8/12/16

%I would like to thank John Donoghue, Basem El-Menoufi, Panayotis Kevrekidis, and William ?Bill? Barnes 
%for useful conversations and inspiration related to this project. This work has been supported in part 
%by the National Science Foundation under grants NSF PHY15-20292 and NSF PHY12-25915.

function slope = slope_btwn(f,next,prev)
%returns slope at point between next and prev

dx = f(2,prev+1)-f(2,prev);

slope = (f(1,next)-f(1,prev))/2/dx;

end

