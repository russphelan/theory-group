%Author: Russell J. Phelan 
%Date: 10-19-15

function [a_prime] = equations(t,a,expandContract,simType,t0)
%implements equations for each simulation type

%global parameters
swError = 'The Expand/Contract switch parameter to the Friedmann equation was set improperly';
p = 1/6/pi/t0^2; %energy density constant. In the case of stationary dust, this has only a mass contribution. 

switch simType
    case 0
%friedman equation governs evolution of scale factor with time. set up in
%case of pressureless dust. 
%requires a>0. Program will terminate if a_prime values become complex. 
g = 1; %Newton's gravitational constant

if expandContract==1
    a_prime = sqrt(8*pi*g*p/a);
elseif expandContract==0
    a_prime = -sqrt(8*pi*g*p/a);
else 
    display(swError);
end
return;
    case 1
%friedman equation with quantum correction. Correction is calculated as
%integral by hand using the classical analytic ODE solution as the
%differential part of the integro-differential equation.

%parameters:
Ns = 1;
Mr = 1;
C1 = log(-1*Mr*t)/t;
C2 = (log(-1*Mr*t)+1)/t^2;


if expandContract==1
    a_prime = sqrt( 1/a*(8*pi*p/3 + Ns/2430/pi*(19*C1/t0^2/t + 26*C2/t0^2)));
elseif expandContract==0
    a_prime = -1*sqrt( 1/a*(8*pi*p/3 + Ns/2430/pi*(19*C1/t0^2/t + 26*C2/t0^2)));
else 
    display(swError);
end
return;
    case 2
        display('The full quantum simulation has not yet been implemented');
        return;
end
end

