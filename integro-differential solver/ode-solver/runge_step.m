%Author: Russell J. Phelan 
%Date: 9-9-16

%I would like to thank John Donoghue, Basem El-Menoufi, Panayotis Kevrekidis, and William ?Bill? Barnes 
%for useful conversations and inspiration related to this project. This work has been supported in part 
%by the National Science Foundation under grants NSF PHY15-20292 and NSF PHY12-25915.

function [next_func_array] = runge_step(func_array,i,area,stepSize,expand_or_contract,simType)
%Takes a function as a 2xn array, returns the function with one more
%Runge-Kutta solver stepSize taken. 

if ~isnan(func_array(1,i))
    next_func_array = func_array;
    return;
end
        
assert((i <= size(func_array,2)),'attempt to overfill function object')


old_t = func_array(2,i-1);
old_a = func_array(1,i-1);
t0 = func_array(2,4);

k1 = equations(old_t,old_a,expand_or_contract,simType,t0,area);
k2 = equations(old_t + stepSize/2,old_a + (stepSize/2)*k1,expand_or_contract,simType,t0,area);
k3 = equations(old_t + stepSize/2,old_a + (stepSize/2)*k2,expand_or_contract,simType,t0,area);
k4 = equations(old_t + stepSize,old_a + stepSize*k3,expand_or_contract,simType,t0,area);
  
if(~isreal(k1) || ~isreal(k2) || ~isreal(k3) || ~isreal(k4))
    %display('things got complex')
    next_a = NaN; %real(old_a + (stepSize/6)*(k1 + 2*k2 + 2*k3 + k4)); 
else
    next_a = old_a + (stepSize/6)*(k1 + 2*k2 + 2*k3 + k4);
end

func_array(1,i) = next_a;
next_func_array = func_array;
end

