%Author: Russell J. Phelan 
%Date: 10-19-15

clear all; 

%user dialog, initial conditions
t0 = input('Enter initial time: (seconds)');
endInt = input('Enter end time: (seconds)');
a0 = input('Enter intial scale factor value: (dimensionless)');
sw = input('Type 1 for expanding universe, 0 for contracting universe.');

step = .01; %step time in seconds

totalSteps = (endInt-t0)/step;

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
   
    k1 = friedman(oldT,old_a,sw);
    k2 = friedman(oldT + step/2,old_a + (step/2)*k1,sw);
    k3 = friedman(oldT + step/2,old_a + (step/2)*k2,sw);
    k4 = friedman(oldT + step,old_a + step*k3,sw);
   
    nextT = oldT + step;
    next_a = old_a + (step/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    if imag(next_a)~=0
        break;
    end
   
    a_array = [a_array next_a];
    t_array = [t_array nextT];
    
    old_a = next_a;
    oldT = nextT;
end

%derivative calculations
aPrime_array = derivOutputArray(a_array,t_array);

aDoublePrime_array = derivOutputArray(aPrime_array,t_array);

%plotting
subplot(2,2,1);

plot(t_array,a_array,'LineWidth',2);

%labels
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');

subplot(2,2,2);
plot(t_array,aPrime_array,'LineWidth',2);

xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\dot{a}(t)$','FontSize',14,'interpreter','latex');
title('Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

subplot(2,2,3);
plot(t_array,aDoublePrime_array,'LineWidth',2);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
title('Second Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');



