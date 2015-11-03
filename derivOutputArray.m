%Author: Russell J. Phelan
%Date: 11-3-15

function y' = derivOutputArray(y,x)
%derivative function takes two arrays to represent function, returns one
%array, the y' prime array. Can be matched with same time array as used for
%original function. 

y' = [];

%difference quotient

dx = x(3)-x(1); %two times steps chosen to find slope symmetric about given point

for i=2:
    
    derivAtPoint = (y(i-1)-y(i+1))/dx

    y'(i) = derivAtPoint;

end

