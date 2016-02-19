function [next_func_array] = runge_step( func_array,stepSize,sw,simType)
%Takes a function as a 2xn array, returns the function with one more
%Runge-Kutta solver stepSize taken. 

i = 1;
while (i<= size(func_array,2) && ~isnan(func_array(1,i)))
    i = i+1;
end

if(i > size(func_array,2))
    disp('Function object is filled')
    return;
end

old_t = func_array(2,i-1);
old_a = func_array(1,i-1);
t0 = func_array(2,1);

k1 = equations(old_t,old_a,sw,simType,t0);
k2 = equations(old_t + stepSize/2,old_a + (stepSize/2)*k1,sw,simType,t0);
k3 = equations(old_t + stepSize/2,old_a + (stepSize/2)*k2,sw,simType,t0);
k4 = equations(old_t + stepSize,old_a + stepSize*k3,sw,simType,t0);
  
next_a = old_a + (stepSize/6)*(k1 + 2*k2 + 2*k3 + k4);
func_array(1,i) = next_a;
next_func_array = func_array;
end

