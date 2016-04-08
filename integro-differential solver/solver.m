%Author: Russell J Phelan
%Date: 4/11/16
%Solver script for integro-differential equations

clear all; 

%THIS SHOULD BE CLEANED UP WHEN YOU'RE DONE WITH EVERYTHING ELSE

%user dialog, initial conditions
simType = input('Choose simulation type: \n(0)Classical Friedmann \n(1)Quantum Approximation \n(2)Full Quantum Behavior');


t0 = input('Enter initial time: (seconds)');
endInt = input('Enter end time: (seconds)');
a0 = input('Enter intial scale factor value: (dimensionless)');
sw = input('Type 1 for expanding universe, 0 for contracting universe.');
eps_steps = 1;
parab_size = 1;
epsilon = .01;
N = 100;

%setting up steps
step = .001; %in seconds
totalSteps = (endInt-t0)/step;
totalSteps = round(totalSteps) + 1;

%setting up scale factor a(t) function object, and derivatives
[scale_factor,scale_1deriv,scale_2deriv,r_func1,r_func2,r_func3,basem_scale_factor] = deal(NaN(2,totalSteps));
for i=1:totalSteps
    [scale_factor(2,i),scale_1deriv(2,i),scale_2deriv(2,i),basem_scale_factor(2,i),r_func1(2,i),r_func2(2,i),r_func3(2,i)] = deal((i-1)*step);
end
scale_factor(1,1) = a0; %initial conditions
basem_scale_factor(1,1) = a0; %initial conditions
scale_factor(2,1) = t0; 
area_matrix = [];

%performing iterations for classical behavior
for curr_t_index=1:totalSteps-1
    
    %calculate next step of R function. going to be passed to
    %causal_nonlocal
%     r_func = r_funcs(r_func,curr_t_index,scale_factor,scale_1deriv,scale_2deriv,1);
%     
%     if(curr_t_index>20)
%         %integrate from t=t0 to current t
%         area = causal_nonlocal_int(1, curr_t_index, r_func, epsilon, step, parab_size);
%         area_matrix = [area_matrix area];
%     end
    
    %calculate next runge-kutta step, update array
    scale_factor = runge_step(scale_factor,0,step,sw,simType); %runge-kutta algorithm
    
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

%performing iterations to reproduce Basem's plots

for curr_t_index=1:totalSteps-1
    
    %calculate first r_func from classical vals
    r_func1 = r_funcs(r_func1,curr_t_index,scale_factor,scale_1deriv,scale_2deriv,1);
    
    %calculate second r_func from classical vals
    r_func2 = r_funcs(r_func2,curr_t_index,scale_factor,scale_1deriv,scale_2deriv,2);
    
    %calculate third r_func from classical vals
    r_func3 = r_funcs(r_func3,curr_t_index,scale_factor,scale_1deriv,scale_2deriv,3);
    
    %calculate area for current t val from the first r_func
    if(curr_t_index>20)
         %integrate from t=t0 to current t for each of the three integrals
         area1 = causal_nonlocal_int(1, curr_t_index, r_func1, epsilon, step, parab_size);
         area2 = causal_nonlocal_int(1, curr_t_index, r_func2, epsilon, step, parab_size);
         area3 = causal_nonlocal_int(1, curr_t_index, r_func3, epsilon, step, parab_size);
         
         %set coefficients for combining areas
         coef1 = 6*sqrt(scale_factor(1,curr_t_index))*scale_2deriv(1,curr_t_index-1);
         coef2 = 6*scale_1deriv(1,curr_t_index-1)/sqrt(scale_factor(1,curr_t_index));
         coef3 = 12*sqrt(scale_factor(1,curr_t_index))*scale_1deriv(1,curr_t_index-1);
         
         %the actual combining of areas
         area = N*(coef1*area1 + coef2*area2 + coef3*area3); 
         
         %storing area for plotting
         area_matrix = [area_matrix area]; 
         
         %calculate next runge_step using the area just calculated
         basem_scale_factor = runge_step(basem_scale_factor,area,step,sw,2); %runge-kutta algorithm
    else
         basem_scale_factor = runge_step(basem_scale_factor,0,step,sw,2); %runge-kutta algorithm
    end
end

%plotting
lw = 1; %sets linewidth for all plots

%scale factor plot
subplot(2,3,1);
plot(scale_factor(2,1:totalSteps),scale_factor(1,1:totalSteps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');

%first derivative plot
subplot(2,3,2);
plot(scale_1deriv(2,1:totalSteps),scale_1deriv(1,1:totalSteps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\dot{a}(t)$','FontSize',14,'interpreter','latex');
title('Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

%second derivative plot
subplot(2,3,3);
plot(scale_2deriv(2,1:totalSteps),scale_2deriv(1,1:totalSteps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
title('Second Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

%area plot
subplot(2,3,4);
plot(scale_2deriv(2,1:length(area_matrix)),area_matrix(1,:),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
title('Area','FontSize',18,'FontWeight','bold','interpreter','latex');

%basem scale factor plot
subplot(2,3,5);
plot(basem_scale_factor(2,1:totalSteps),basem_scale_factor(1,1:totalSteps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Basem Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');
