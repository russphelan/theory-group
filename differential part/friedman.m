%Author: Russell J. Phelan 
%Date: 10-19-15

function [a_prime] = friedman(t,a,sw)
%friedman equation governs evolution of scale factor with time. set up in
%case of pressureless dust. 

%requires a>0. 

%This is a quadratic ODE. 

p = 1; %energy density constant. In the case of stationary dust, this has only a mass contribution. 
g = 1; %Newton's gravitational constant

if sw==1
    a_prime = sqrt(8*pi*g*p/a);
elseif sw==0
    a_prime = -sqrt(8*pi*g*p/a);
else 
    display('The switch parameter to the Friedmann equation was set improperly')
end
end

