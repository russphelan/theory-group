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
eps_steps = 1;
parab_size = 1;
epsilon = .01;

%setting up steps
step = .001; %in seconds
totalSteps = (endInt-t0)/step;
totalSteps = round(totalSteps) + 1;

%setting up scale factor a(t) function object, and derivatives
[scale_factor,scale_1deriv,scale_2deriv,r_func] = deal(NaN(2,totalSteps));
for i=1:totalSteps
    [scale_factor(2,i),scale_1deriv(2,i),scale_2deriv(2,i),r_func(2,i)] = deal((i-1)*step);
end
scale_factor(1,1) = a0; %initial conditions
scale_factor(2,1) = t0; 
area_matrix = [];

%performing iterations
for curr_t_index=1:totalSteps-1
    
    %calculate next step of R function. going to be passed to
    %causal_nonlocal
    r_func = r_funcs(r_func,curr_t_index,scale_factor,scale_1deriv,scale_2deriv,1);
    
    if(curr_t_index>20)
        %integrate from t=t0 to current t
        area = causal_nonlocal_int(1, curr_t_index, r_func, epsilon, step, parab_size);
        area_matrix = [area_matrix area];
    end
    
    %calculate next runge-kutta step, update array
    scale_factor = runge_step(scale_factor,step,sw,simType); %runge-kutta algorithm
    
    %calculate next set of derivative steps, update arrays
    if(curr_t_index>1)
        scale_1deriv(1,curr_t_index) = slope_btwn(scale_factor,curr_t_index,curr_t_index-1);
    end
    if(curr_t_index>3)
        scale_2deriv(1,curr_t_index) = slope_btwn(scale_1deriv,curr_t_index-1,curr_t_index-2);
    end
    
    %error checks: 
    assert(scale_factor(1,curr_t_index)>=0, 'Behavior is not well-defined for a negative scale factor');
    assert(~isnan(scale_factor(1,curr_t_index)),'The scale factor is not a number');
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

%second derivative plot
subplot(2,2,4);
plot(scale_2deriv(2,1:length(area_matrix)),area_matrix(1,:),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
title('Area','FontSize',18,'FontWeight','bold','interpreter','latex');