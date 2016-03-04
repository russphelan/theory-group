%Author: Russell J Phelan
%Date: 12/8/15
%Solver script for integro-differential equations

clear all; 

%global settings
step = .001; %in seconds


%THIS SHOULD BE CLEANED UP WHEN YOU'RE DONE WITH EVERYTHING ELSE

%user dialog, initial conditions
simType = input('Choose simulation type: \n(0)Classical Friedmann \n(1)Quantum Approximation \n(2)Full Quantum Behavior');
assert(simType ~= 2,'Full quantum correction not yet implemented');

t0 = input('Enter initial time: (seconds)');
endInt = input('Enter end time: (seconds)');
a0 = input('Enter intial scale factor value: (dimensionless)');
sw = input('Type 1 for expanding universe, 0 for contracting universe.');

%setting up steps
totalSteps = (endInt-t0)/step;
totalSteps = round(totalSteps) + 1; %+1 ensures evaluation at end time given...I know it looks strange.

%setting up scale factor a(t) function object, and derivatives
[scale_factor,scale_1deriv,scale_2deriv] = deal(NaN(2,totalSteps)); %fill them with NaNs
for i=1:totalSteps
    [scale_factor(2,i),scale_1deriv(2,i),scale_2deriv(2,i)] = deal((i-1)*step); %fill second row with t values spaced with 'step'
end
scale_factor(1,1) = a0; %initial conditions
scale_factor(2,1) = t0; 

%performing iterations
for i=1:totalSteps-1
    
    %calculate next runge-kutta step, update array
    scale_factor = runge_step(scale_factor,step,sw,simType); %runge-kutta algorithm
    
    %calculate next set of derivative steps, update arrays
    if(i>1)%filled with NaNs until i=2
        scale_1deriv(1,i) = slope_btwn(scale_factor,i,i-1);
    end
    if(i>3) %filled with NaNs till i=4
        scale_2deriv(1,i) = slope_btwn(scale_1deriv,i-1,i-2);
    end
    
    %error checks: 
    assert(scale_factor(1,i)>=0, 'Behavior is not well-defined for a negative scale factor');
    assert(~isnan(scale_factor(1,i)),'The scale factor is not a number');
end

%plotting
lw = 1; %sets linewidth for all plots

%scale factor plot
subplot(2,2,1);
plot(scale_factor(2,1:totalSteps),scale_factor(1,1:totalSteps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');

%first derivative plot
subplot(2,2,2);
plot(scale_1deriv(2,1:totalSteps),scale_1deriv(1,1:totalSteps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\dot{a}(t)$','FontSize',14,'interpreter','latex');
title('Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

%second derivative plot
subplot(2,2,3);
plot(scale_2deriv(2,1:totalSteps),scale_2deriv(1,1:totalSteps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
title('Second Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');