%Author: Russell J. Phelan 
%Date: 10-5-15

function area = causal_nonlocal_int(startInt,endInt,e,t,dx)
%Causal non-local function implemented with small epsilon in place of
%limiting process. dx must be smaller than e. e<<1 is also required. 

%simpson's rule implementation. here, dx is what we call dtPrime in the
%paper.

%setting bounds on integral
if endInt<=t-e
    numInts = (endInt-startInt)/dx;
else
    numInts = ((t-e)-startInt)/dx;
end

area = 0;

for i=0:(numInts-1)
    a = startInt + dx*i;
    b = startInt + dx*i + dx;
    
    area = area + ((b-a)/6)*(f(a)/(t-a) + 4*f((a+b)/2)/(t-(a+b)/2) + f(b)/(t-b)); 
end

%implements delta function term to compensate for divergence
area = area + log(e)*f(t);

end

