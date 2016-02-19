%Author: Russell J. Phelan 
%Date: 10-5-15

function area = causal_nonlocal_int(startInt,endInt,scale_factor,endIndex,scale_1deriv,scale_2deriv,e,t,step,widthSlice)
%Causal non-local function implemented with small epsilon in place of
%limiting process. prints error. dx must be smaller than e. e<<1 is also required. 

%INTEGRAL CALCULATION
%simpson's rule implementation. here, dx is what we call dtPrime in the
%paper.

dx = step*widthSlice;
%setting bounds on integral
if endInt<=t-e
    numInts = (endInt-startInt)/dx;
else
    numInts = ((t-e)-startInt)/dx;
end

numInts = round(numInts);
area = 0;

for i=0:1:(numInts-1)
    %a and b are current values of tPrime. They mark the edges of the
    %current slice of integral. 
    a = startInt + dx*i;
    b = startInt + dx*i + dx;
    
    aIndex = round((a-startInt)/step + 1);
    bIndex = round((b-startInt)/step + 1);
    halfIndex = round(((a+b)/2 - startInt)/step + 1);
    
   area = area + ((b-a)/6)*(f(aIndex, scale_factor, scale_1deriv, scale_2deriv, 1)/(t-a) + 4*f(halfIndex, scale_factor, scale_1deriv, scale_2deriv, 1)/(t-(a+b)/2) + f(bIndex, scale_factor, scale_1deriv, scale_2deriv, 1)/(t-b)); 
end


%implements delta function term to compensate for divergence
tIndex = round((t-startInt)/step + 1);
area = area + log(e)*f(tIndex, scale_factor, scale_1deriv, scale_2deriv, 1); 
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