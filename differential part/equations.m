%Author: Russell J. Phelan 
%Date: 10-19-15

function [a_prime] = equations(t,a,expandContract,simType,t0)

%error conditions
swError = 'The Expand/Contract switch parameter to the Friedmann equation was set improperly. Needs 1 or 0';
zeroError = 'The quantum correction behavior is not defined for t=0';
 
assert((expandContract == 0) || (expandContract == 1),swError);

%implements equations for each simulation type
switch simType
    
    case 0
        
%CLASSICAL BEHAVIOR
%requires a>0. Program will terminate if a_prime values become complex. 
g = 1; %Newton's gravitational constant
p = 1; %density of dust

if expandContract==1
    a_prime = sqrt(8*pi*g*p/a);
elseif expandContract==0
    a_prime = -sqrt(8*pi*g*p/a);
end
return;

    case 1
        
%QUANTUM CORRECTION, HAND INTERGRALS
%error checks
assert(t0 ~= 0,zeroError);

%parameters:
p = 1/6/pi/t0^2; %energy density
Ns = 1;
Mr = 1;
C1 = log(Mr*abs(t))/t;
C2 = (log(Mr*abs(t))+1)/t^2;

if expandContract==1
    assert(false, 'expansion behavior is not implemented for the quantum correction');
elseif expandContract==0
    a_prime = -1*sqrt( 1/a*(8*pi*p/3 + Ns/2430/pi*(19*C1/t0^2/t + 26*C2/t0^2)));
    if ~isreal(a_prime)
        disp(sprintf('the value for :a: is: %d /n',a))
        disp(sprintf('the value for :gamma: is: %d /n',(8*pi*p/3 + Ns/2430/pi*(19*C1/t0^2/t + 26*C2/t0^2))))
        disp(sprintf('the value for :alpha: is: %d /n',8*pi*p/3))
        disp(sprintf('the value for :beta: is: %d /n',Ns/2430/pi*(19*C1/t0^2/t + 26*C2/t0^2)))
        disp(sprintf('the value for :firstlog: is: %d /n',19*C1/t0^2/t))
        disp(sprintf('the value for :secondlog: is: %d /n',26*C2/t0^2))
        assert(isreal(a_prime),'The scale factor deriv has become imaginary');
    end
end
return;

    case 2
        %NOT YET IMPLEMENTED
end
end

