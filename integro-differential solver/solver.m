%Author: Russell J Phelan
%Date: 12/8/15
%Solver script for integro-differential equations

clear all; 

%THIS SHOULD BE CLEANED UP WHEN YOU'RE DONE WITH EVERYTHING ELSE

%user dialog, initial conditions
simType = input('Choose simulation type: \n(0)Classical Friedmann \n(1)Quantum Approximation \n(2)Full Quantum Behavior');

assert(simType ~= 2,'Full quantum correction not yet implemented');

t0 = input('Enter initial time: (seconds)');
endInt = input('Enter end time: (seconds)');
a0 = input('Enter intial scale factor value: (dimensionless)');
sw = input('Type 1 for expanding universe, 0 for contracting universe.');

%setting up steps
step = .001; %in seconds
totalSteps = (endInt-t0)/step;
totalSteps = round(totalSteps) + 1;

%setting up function object
scale_factor = NaN(2,totalSteps);
for i=1:totalSteps
    scale_factor(2,i) = (i-1)*step;
end
scale_factor(1,1) = a0; %initial conditions
scale_factor(2,1) = t0; 

%performing iterations
for i=1:totalSteps-1
    
    scale_factor = runge_step(scale_factor,step,sw,simType); %runge-kutta algorithm
    
    scale_factor(1,2);
    %error checks: 
    assert(scale_factor(1,i)>=0, 'Behavior is not well-defined for a negative scale factor');
    assert(~isnan(scale_factor(1,i)),'The scale factor is not a number');
end



% %DERIVATIVE SECTION NEEDS TO BE REWRITTEN USING SLOPE_BTWN FUNCTION
% 
% %derivative calculations
% aPrime_array = derivOutputArray(a_array,t_array);
% aDoublePrime_array = derivOutputArray(aPrime_array,t_array);



%plotting

lw = 1; %sets linewidth for all plots
subplot(2,2,1);

plot(scale_factor(2,1:totalSteps),scale_factor(1,1:totalSteps),'LineWidth',lw);
 
%labels
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');



%DERIVATIVE PLOTS NOT FOR USE YET
% 
% subplot(2,2,2);
% plot(t_array,aPrime_array,'LineWidth',lw);
% 
% xlabel('Time (s)','FontSize',14,'interpreter','latex');
% ylabel('$\dot{a}(t)$','FontSize',14,'interpreter','latex');
% title('Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');
% 
% subplot(2,2,3);
% plot(t_array,aDoublePrime_array,'LineWidth',lw);
% xlabel('Time (s)','FontSize',14,'interpreter','latex');
% ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
% title('Second Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');