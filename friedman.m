%Author: Russell J. Phelan 
%Date: 10-19-15

function [a_prime] = friedman(t,a)
%friedman equation governs evolution of scale factor with time. set up in
%case of pressureless dust. 

%This is a quadratic ODE, I am taking the positive root for now. Reduction of order could be appropriate.

p = 1; %energy density function. In the case of stationary dust, this has only a mass contribution. 
g = 1; %Newton's gravitational constant

a_prime = sqrt(8*pi*g*p/a);
end

