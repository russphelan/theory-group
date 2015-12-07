%Author: Russell J. Phelan 
%Date: 10-19-15

clear all; 

%user dialog, initial conditions
simType = input('Choose simulation type: \n(0)Classical Friedmann \n(1)Quantum Approximation \n(2)Full Quantum Behavior');
assert(simType ~= 2,'Full quantum correction not yet implemented');
t0 = input('Enter initial time: (seconds)');
endInt = input('Enter end time: (seconds)');
a0 = input('Enter intial scale factor value: (dimensionless)');
sw = input('Type 1 for expanding universe, 0 for contracting universe.');

step = .001; %step size in seconds
totalSteps = (endInt-t0)/step;
totalSteps = round(totalSteps);

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
    
    %error checks: 
    assert(old_a>=0, 'behavior is not defined for a negative scale factor');
    assert(~isnan(next_a),'The scale factor is not a number');
    
    k1 = equations(oldT,old_a,sw,simType,t0);
    k2 = equations(oldT + step/2,old_a + (step/2)*k1,sw,simType,t0);
    k3 = equations(oldT + step/2,old_a + (step/2)*k2,sw,simType,t0);
    k4 = equations(oldT + step,old_a + step*k3,sw,simType,t0);
   
    nextT = oldT + step;
    next_a = old_a + (step/6)*(k1 + 2*k2 + 2*k3 + k4);
   
    a_array = [a_array next_a];
    t_array = [t_array nextT];
    
    old_a = next_a;
    oldT = nextT;
end

%derivative calculations
aPrime_array = derivOutputArray(a_array,t_array);
aDoublePrime_array = derivOutputArray(aPrime_array,t_array);

%plotting

lw = 1; %sets linewidth for all plots
subplot(2,2,1);

plot(t_array,a_array,'LineWidth',lw);

%labels
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');

subplot(2,2,2);
plot(t_array,aPrime_array,'LineWidth',lw);

xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\dot{a}(t)$','FontSize',14,'interpreter','latex');
title('Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

subplot(2,2,3);
plot(t_array,aDoublePrime_array,'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
title('Second Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');