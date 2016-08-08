%Author: Russell J Phelan
%Date: 6/12/16
%Solver script for integro-differential equations

clear all; 

%INITIAL CONDITIONS/PARAMETERS
%simType = input('Choose simulation type: \n(0)Classical Friedmann \n(1)Quantum Approximation \n(2)Full Quantum Behavior');
simType = 0; %for testing
t0 = input('Enter initial time: (seconds)');
tf = input('Enter end time: (seconds)');
a0 = input('Enter intial scale factor value: (dimensionless)');
%a0 = 1; %for testing
%expand_or_contract = input('Type 1 for expanding universe, 0 for contracting universe.'); %chooses which branch of scale factor equations to use
expand_or_contract = 0; %for testing

rect_thickness = .001;
e = .001;
N = 100; %number of scalars, also N_s in paper. 
step = .001; %in seconds
total_steps = (tf-t0)/step;
total_steps = ceil(total_steps) + 1;

%setting up functions filled with NaNs, and t values, spaced by the step. 
%NaNs are first row, t values are second row. 
[scale_factor,scale_1deriv,scale_2deriv,scale_3deriv,r_func1,r_func2,r_func3,basem_scale_factor,basem_scale_1deriv] = deal(NaN(2,total_steps));
for i=1:total_steps
    [scale_factor(2,i),scale_1deriv(2,i),scale_2deriv(2,i),scale_3deriv(2,i),basem_scale_factor(2,i),basem_scale_1deriv(2,i),r_func1(2,i),r_func2(2,i),r_func3(2,i)] = deal(t0 + (i-1)*step);
end

scale_factor(1,1) = a0; %initial conditions
basem_scale_factor(1,1) = a0; %initial conditions
area_matrix = [0 0 0];

%CLASSICAL ITERATIONS
for curr_t_index=1:total_steps
    
    %calculate next runge-kutta step, update array, 0 is for no area
    scale_factor = runge_step(scale_factor,curr_t_index,0,step,expand_or_contract,simType); %runge-kutta algorithm
    
    %calculate next set of derivative steps, update arrays
    if(curr_t_index>=3)
        scale_1deriv(1,curr_t_index) = slope_btwn(scale_factor,curr_t_index,curr_t_index-2);
    end
    if(curr_t_index>=5)
        scale_2deriv(1,curr_t_index) = slope_btwn(scale_1deriv,curr_t_index,curr_t_index-2);
    end
    if(curr_t_index>=7)
        scale_3deriv(1,curr_t_index) = slope_btwn(scale_2deriv,curr_t_index,curr_t_index-2);
    end
    
    %ERROR CHECKS
    assert(scale_factor(1,curr_t_index)>=0, 'Behavior is not well-defined for a negative scale factor');
    
    %stops iterations when scale factor becomes NaN, but you can still see
    %the plot thus far. Notifies user. 
    if(isnan(scale_factor(1,curr_t_index)))
        display('scale factor is NaN at last iteration')
        break;
    end
    
    %stops iterations when scale factor becomes complex, letting you see
    %the plots thus far. Notifies user. 
    if(~isreal(scale_factor(1,curr_t_index)))
        display('scale factor is complex at last iteration')
        break;
    end
end

%NUMERICAL AREA
for t=1:total_steps
    
    
    
    %calculate area up to current t value
    if(t>=7)
        %r functions
        r_func1 = r_funcs(r_func1,t,scale_factor,scale_1deriv,scale_2deriv,scale_3deriv,1);
        r_func2 = r_funcs(r_func2,t,scale_factor,scale_1deriv,scale_2deriv,scale_3deriv,2);
        r_func3 = r_funcs(r_func3,t,scale_factor,scale_1deriv,scale_2deriv,scale_3deriv,3);
        
        %integrate from t=t0 to current t for each of the three integrals
        area1 = causal_nonlocal_int(t, r_func1, e, step, rect_thickness);
        area2 = causal_nonlocal_int(t, r_func2, e, step, rect_thickness);
        area3 = causal_nonlocal_int(t, r_func3, e, step, rect_thickness);
        
        %set coefficients for combining areas
        coef1 = 6*sqrt(scale_factor(1,t))*scale_2deriv(1,t);
        coef2 = 6*scale_1deriv(1,t)^2/sqrt(scale_factor(1,t));
        coef3 = 12*sqrt(scale_factor(1,t))*scale_1deriv(1,t);
        
        %the actual combining of areas
        area = N*(coef1*area1 + coef2*area2 + coef3*area3);
        
        %storing area for plotting
        area_matrix = [area_matrix area];
        
        %calculate next runge_step using the area just calculated
        basem_scale_factor = runge_step(basem_scale_factor,t,area,step,expand_or_contract,2); %runge-kutta algorithm
        basem_scale_1deriv(1,t) = slope_btwn(basem_scale_factor,t,t-2);
        
    else
        basem_scale_factor = runge_step(basem_scale_factor,t,0,step,expand_or_contract,2); %runge-kutta algorithm
        if t>=3
            basem_scale_1deriv(1,t) = slope_btwn(basem_scale_factor,t,t-2);
        end
    end
end

%Basem's hand-calculated area function (for comparison) 

%Expanding phase function
%Parameters
N_s = 1;
Mu_r = 1; 

%t array initialization
for i=1:total_steps
    t(i) = deal(t0 + (i-1)*step);
end

%basem_area = -N_s/2430/pi*(19*(log(Mu_r*t)+log(t/t0-1))./t.^2/t0^2+26*(log(Mu_r*t)+log(t/t0-1)+(t/t0-1))./t.^2/t0^2);

%N_s/2430/pi*(19*(log(Mu_r*a)+log(a/t0-1))/a^2/t0^2+26*(log(Mu_r*a)+log(a/t0-1)+(a/t0-1))/a^2/t0^2)

%Contracting phase function
%basem_area = -N_s/2430/pi*(19*log(-Mu_r.*t)/t0^2./t.^2 + 26*(log(-Mu_r.*t)+1)/t0^2./t.^2);

%plotting
lw = 1; %sets linewidth for all plots

%scale factor plot
 subplot(2,3,1);
 plot(scale_factor(2,1:total_steps),scale_factor(1,1:total_steps),'LineWidth',lw);
 xlabel('Time (s)','FontSize',14,'interpreter','latex');
 ylabel('$a(t)$','FontSize',14,'interpreter','latex');
 title('Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');

%first derivative plot
subplot(2,3,2);
plot(scale_1deriv(2,1:total_steps),scale_1deriv(1,1:total_steps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$\dot{a}(t)$','FontSize',14,'interpreter','latex');
title('Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

%second derivative plot
% subplot(2,3,3);
% plot(scale_2deriv(2,1:total_steps),scale_2deriv(1,1:total_steps),'LineWidth',lw);
% xlabel('Time (s)','FontSize',14,'interpreter','latex');
% ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
% title('Second Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

%third derivative plot
% subplot(2,3,4);
% plot(scale_3deriv(2,1:total_steps),scale_3deriv(1,1:total_steps),'LineWidth',lw);
% xlabel('Time (s)','FontSize',14,'interpreter','latex');
% ylabel('$\ddot{a}(t)$','FontSize',14,'interpreter','latex');
% title('Second Derivative','FontSize',18,'FontWeight','bold','interpreter','latex');

%area plot
subplot(2,3,3);
plot(scale_2deriv(2,1:length(area_matrix)),area_matrix(1,:),'LineWidth',lw,'Color','r');
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('Area from $t_0$','FontSize',14,'interpreter','latex');
title('Area','FontSize',18,'FontWeight','bold','interpreter','latex');

%basem scale factor plot
subplot(2,3,4);
plot(basem_scale_factor(2,1:total_steps),basem_scale_factor(1,1:total_steps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Basem Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');

%basem scale factor derivative plot
subplot(2,3,5);
plot(basem_scale_1deriv(2,1:total_steps),basem_scale_1deriv(1,1:total_steps),'LineWidth',lw);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('$a(t)$','FontSize',14,'interpreter','latex');
title('Basem Scale Factor','FontSize',18,'FontWeight','bold','interpreter','latex');

%basem area comparison function, contraction
% subplot(2,3,6);
% plot(t,basem_area,'LineWidth',lw);
% xlabel('Time (s)','FontSize',14,'interpreter','latex');
% ylabel('Area from $-\infty$','FontSize',14,'interpreter','latex');
% title('Basem Area ','FontSize',18,'FontWeight','bold','interpreter','latex');