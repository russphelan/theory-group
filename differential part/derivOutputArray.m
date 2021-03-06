%Author: Russell J. Phelan
%Date: 11-3-15

function yPrime = derivOutputArray(y,x)
%Arugments: y is array of y(x) values. x is array of x values. They must be
%the same size. 

%Returns: y' array of y'(x) values. Matches with original x values. 

%error handling for input arguments
assert(length(y)==length(x),'the argument arrays are not of equal length');

%init output array
yPrime = zeros([1 length(x)]);

%difference quotient

dx = x(3)-x(1); %two time steps chosen to find slope symmetric about given point

for i=2:length(x)-1
    derivAtPoint = (y(i+1)-y(i-1))/dx;

    yPrime(i) = derivAtPoint;
end

%fix beginning and end of output array
yPrime(1) = yPrime(2);
yPrime(length(x)) = yPrime(length(x)-1);

end

