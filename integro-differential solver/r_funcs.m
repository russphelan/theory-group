%Author: Russell J. Phelan 
%Date: 10-5-15

function new_r_func = r_funcs(r_func,curr_t_index,scale_factor,scale_1deriv,scale_2deriv,which_r)
%function for multiplying by script L. Stores all of the R functions. 


alpha=1;
beta=1;
gamma=1;

deriv_index = curr_t_index-1;


if curr_t_index<=4
    r_func(1,curr_t_index) = 0;
else
    switch(which_r)
        case 1 %for R1
            
        r_func(1,curr_t_index) = -1*sqrt(scale_factor(1,curr_t_index))*scale_2deriv(1,deriv_index)*(6*alpha + 2*beta + 2*gamma) ...
            - scale_1deriv(1,deriv_index)^2/sqrt(scale_factor(1,curr_t_index))*(6*alpha + beta);
        
        assert(~isnan(r_func(1,curr_t_index)),'this r_funcs value is NaN')
       
        case 2 %for R2
            
            r_func(1,curr_t_index) = -1*sqrt(scale_factor(1,curr_t_index))*scale_2deriv(1,deriv_index)*(12*alpha + 1*beta + -2*gamma) ...
            - scale_1deriv(1,deriv_index)^2/sqrt(scale_factor(1,curr_t_index))*(12*alpha + 5*beta + 6*gamma);
       
        case 3 %for R3  
            
            r_func(1,curr_t_index) = sqrt(scale_factor(1,curr_t_index))*scale_2deriv(1,deriv_index)*(6*alpha + 2*beta + 2*gamma) ...
            + scale_1deriv(1,deriv_index)^2/sqrt(scale_factor(1,curr_t_index))*(6*alpha + beta);
        
        otherwise
            error('The switch parameter for choosing the R function does not correspond to any of the current functions');
    end
end
new_r_func = r_func;
end

