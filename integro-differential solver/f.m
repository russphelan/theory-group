%Author: Russell J. Phelan 
%Date: 10-5-15

function y = f(index,scale_factor, scale_1deriv, scale_2deriv, whichR)
%function for multiplying by script L. 

%for R1

firstIndex = index-1;
secondIndex = index-2; 

if ((firstIndex==0) || (secondIndex==0))
    y = 0;
    return;
end

alp = 1; %parameters in R1
bet = 1;
gam = 1;

y = -1*sqrt(scale_factor(1,index))*scale_2deriv(1,secondIndex)*(6*alp + 2*bet + 2*gam) - scale_1deriv(1,firstIndex)^2/sqrt(scale_factor(1,index))*(6*alp + bet);
assert(~isnan(y),'y is a NaN');
end

