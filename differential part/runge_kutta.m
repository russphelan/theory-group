%Author: Russell J. Phelan 
%Date: 10-19-15

clear all; 

%user dialog, initial conditions
t0 = input('Enter initial time: (seconds)');
a0 = input('Enter intial scale factor value: (dimensionless)');

step = .01; %step time in seconds

startInt = 0; %start of time interval
endInt = 10;  %end of time interval
totalSteps = (endInt-startInt)/step;

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
t_array = t0;

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

%derivative calculation
aPrime_array = derivOutputArray(a_array,t_array);

%plotting
subplot(1,2,1);

plot(t_array,a_array);

%labels
xlabel('Time (s)','FontSize',12);
ylabel('a(t)','FontSize',12);
title('Scale Factor vs. Time','FontSize',14,'FontWeight','bold');

subplot(1,2,2);
plot(t_array,aPrime_array);

xlabel('Time (s)','FontSize',12);
ylabel('aPrime(t)','FontSize',12);
title('Scale Factor Derivative vs. Time','FontSize',14,'FontWeight','bold');




