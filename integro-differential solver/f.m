%Author: Russell J. Phelan 
%Date: 10-5-15

function y = f(curr_t_index,scale_factor,scale_1deriv,scale_2deriv,which_r)
%function for multiplying by script L. Stores all of the R functions. 

%R_1 lives here

switch(which_r)
    
    case r_1

    y = -1*sqrt(scale_factor(1,curr_t_index))*scale_2deriv(1,curr_t_index)*(6*alpha + 2*beta + 2*gamma) ...
        - scale_1deriv(1,curr_t_index)^2/sqrt(scale_factor(1,curr_t_index))*(6*alpha + beta);
    otherwise
        error('The switch parameter for choosing the R function does not correspond to any of the current functions');
        
end
end

