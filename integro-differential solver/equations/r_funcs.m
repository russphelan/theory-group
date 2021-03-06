%Author: Russell J. Phelan 
%Date: 10-5-15

%I would like to thank John Donoghue, Basem El-Menoufi, Panayotis Kevrekidis, and William ?Bill? Barnes 
%for useful conversations and inspiration related to this project. This work has been supported in part 
%by the National Science Foundation under grants NSF PHY15-20292 and NSF PHY12-25915.

function new_r_func = r_funcs(r_func,t,scale_factor,scale_1deriv,scale_2deriv, scale_3deriv,which_r)
%function for multiplying by script L. Stores all of the R functions. 

alpha=5/11520/pi^2;
beta=-2/11520/pi^2;
gamma=2/11520/pi^2;

switch(which_r)
    case 1 %for R1
        
        r_func(1,t) = -1*sqrt(scale_factor(1,t))*scale_2deriv(1,t)*(6*alpha + 2*beta + 2*gamma) ...
            - scale_1deriv(1,t)^2/sqrt(scale_factor(1,t))*(6*alpha + beta);
        
    case 2 %for R2
        
        r_func(1,t) = -1*sqrt(scale_factor(1,t))*scale_2deriv(1,t)*(12*alpha + beta - 2*gamma) ...
            - scale_1deriv(1,t)^2/sqrt(scale_factor(1,t))*(12*alpha + 5*beta + 6*gamma);
        
    case 3 %for R3
        
        r_func(1,t) = (6*alpha + 2*beta + 2*gamma)*(1/2/sqrt(scale_factor(1,t))*scale_1deriv(1,t)*scale_2deriv(1,t) ...
            + sqrt(scale_factor(1,t))*scale_3deriv(1,t)) + (6*alpha + beta)*(2*scale_1deriv(1,t)*scale_2deriv(1,t)/sqrt(scale_factor(1,t)) ...
            - scale_1deriv(1,t)^2/2/scale_factor(1,t)^(3/2));
        
    case 4 %for testing integrator
        r_func(1,t) = 1/(r_func(2,t)^2);  %cos(r_func(2,t));
        
    otherwise
        error('The switch parameter for choosing the R function does not correspond to any of the current functions');
end

%error conditions
assert(~isnan(r_func(1,t)),'this r_funcs value is NaN at index %u',t)
assert(isreal(r_func(1,t)),'this r_funcs value is complex at index %u',t)
        
new_r_func = r_func;
end