%Author: Russell J. Phelan 
%Date: 6-16-16

function area = causal_nonlocal_int(t, r_func, e, step, rect_thickness)
%gives 0 if not enough space for single rectangle, otherwise returns
%trapazoid rule area with epsilon behavior. 

%Causal non-local function implemented with small epsilon in place of
%limiting process. dx must be smaller than e. e<<1 is also required. 
%t is index where we want to end integration, the current t value 

%step size is integer which determines over how many steps a single
%parabola will approximate

%Causal non-local function implemented with small epsilon in place of
%limiting process. step must be smaller than e. e<<1 is also required. 

%t is an index, gives current t value in calling script when indexed into
%the independent var of a function object

Mu_r = 1;                               %normalization constant
rect_index = ceil(rect_thickness/step); %index interval that corresponds to the thickness of a rectangle, panel, trapezoid, etc when indexed into the independent var of a function object
e_index = ceil(e/step);                 %index interval that corresponds to size of e when used to index into the independent variable of a function object
area = 0;                               %area under the curve. this is returned when the function returns. 
sum = 0;                                %stores a running sum of the integral. does not repesent area under curve until end points are added, and coefficient is multiplied
N = 0;                                  %represents number of panels, or trapezoids used to integrate

%don't start calculating area until t index is late enough to give us room to
%fit in at least one rectangle of the specified thickness
if t < 4+rect_index+e_index
    
    area = 0;
    return;
end

for i=4+1:rect_index:t-e_index-1
    %display(1/(r_func(2,t)-r_func(2,i)))
 
    sum = sum + 2*1/(r_func(2,t)-r_func(2,i))*r_func(1,i); %summing over all terms with 2 coefficient, multiplying each by 1/(t-t')
    N = N+1;
    
end

N = N+1; %would be N+2, since we are about to add two more points, but the number of panels is the number of points-1. Can be understood with a drawing. the points bound the rectangle in terms of the independent variable. 

sum = sum + 1/(r_func(2,t)-r_func(2,4))*r_func(1,4) + 1/(r_func(2,t)-r_func(2,t-e_index))*r_func(1,t-e_index); %adding the first term, and last term, which have no 2 coefficient, multiplying by 1/(t-t')

area = (r_func(2,t)-e-r_func(2,4))/2/N*sum; %multiplying by trapezoid rule coefficient

area = area + log(Mu_r*e)*r_func(1,t); %implements delta function term to compensate for divergence

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