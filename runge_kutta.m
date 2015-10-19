%Author: Russell J. Phelan 
%Date: 10-19-15

clear all; 

step = .001; %step time in seconds

startInt = 0; %start of time interval
endInt = 10;  %end of time interval
totalSteps = (endInt-startInt)/step;

t0 = 0; %initial time value. 
a0 = 1; %initial scale factor

%setting initial conditions
old_a = a0;
oldT = t0;
next_a = 0;
nextT = 0;
k1 = 0;
k2 = 0;
k3 = 0;
k4 = 0;

a_array = a0;
t_array= t0;

%performing iterations
for i=1:totalSteps
   
    k1 = friedman(oldT,old_a);
    k2 = friedman(oldT + step/2,old_a + (step/2)*k1);
    k3 = friedman(oldT + step/2,old_a + (step/2)*k2);
    k4 = friedman(oldT + step,old_a + step*k3);
   
    nextT = oldT + step;
    next_a = old_a + (step/6)*(k1 + 2*k2 + 2*k3 + k4);
   
    a_array = [a_array next_a];
    t_array = [t_array nextT];
    
    old_a = next_a;
    oldT = nextT;
end

%plotting

plot(t_array,a_array);

