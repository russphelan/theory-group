%Author: Russell J. Phelan 
%Date: 10-5-15

function area = causal_nonlocal_int(start_index, end_index, r_func, epsilon,step_size,parab_size)
%Causal non-local function implemented with small epsilon in place of
%limiting process. prints error. dx must be smaller than e. e<<1 is also required. 
%start_index is the index of a function object where we want to start
%integrating the function. end_index works same way. 

%eps_steps is number of steps before end_index that we will integrate to.
%The size of epsilon can be calculated by eps_steps*step where step is the
%step size from "solver.m". 

%step size is integer which determines over how many steps a single
%parabola will approximate

%Causal non-local function implemented with small epsilon in place of
%limiting process. prints error. dx must be smaller than e. e<<1 is also required. 


%INTEGRAL CALCULATION
%simpson's rule implementation. here, dx is what we call dtPrime in the
%paper.

t=r_func(2,end_index); %current t passed to this function


area = 0;

%first integral, with R1
for i=start_index:end_index-1
    a = i;
    mid = i + parab_size;
    b = i + 2*parab_size;
    
    area = area + r_func(1,a)*(r_func(2,b)-r_func(2,a));
    
    %area = area + ((r_func(2,b)-r_func(2,a))/6)*r_func(1,a)/(t-r_func(2,a)) + 4*r_func(1,mid)/(t-r_func(2,mid)) + r_func(1,b)/(t-r_func(2,b))
end

%implements delta function term to compensate for divergence
area = area + log(epsilon)*r_func(1,end_index);


%sprintf('end index is : %f',end_index)
%sprintf('r_func is : %f',r_func(1,end_index))
%assert(~isnan(area),'area went NaN')


% %ERROR CALCULATION
% %searching for a max in the fourth derivative
% small = .01;
% stop = ((t-e)-startInt)/.01;
% stop = round(stop);
% fourthDerivs = [];
% for i=0:1:stop
%     fourthDerivs = [fourth_deriv(startInt+i*small) fourthDerivs];
% end
% fourth_deriv_upperbound = max(fourthDerivs);
% 
% %calculating error based on formula
% error = -(dx^4)/180*((t-e)-startInt)*fourth_deriv_upperbound;
% 
% out = sprintf('The error in calculating the integral is %e',error);
% 
% disp(out)