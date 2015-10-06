%test script for derivative functions

clear all; 


x = 0:.001:2*3.14159;

first = deriv_at_point(x);

y = f(x);

second = second_deriv(x);

third = third_deriv(x);

fourth = fourth_deriv(x);



plot(x,y,'b')
hold on;
plot(x,first,'r')
hold on;
plot(x,second,'g')
hold on;
plot(x,third,'c')
hold on;
plot(x,fourth,'k')

