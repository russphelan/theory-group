%Author: Russell J. Phelan 
%Date: 9-9-16

%I would like to thank John Donoghue, Basem El-Menoufi, Panayotis Kevrekidis, and William ?Bill? Barnes 
%for useful conversations and inspiration related to this project. This work has been supported in part 
%by the National Science Foundation under grants NSF PHY15-20292 and NSF PHY12-25915.

%sets up r_func 
r_funct = deal(NaN(2,total_steps));

for i=1:total_steps
    r_funct(2,i) = deal(t0 + (i-1)*step);
end

for t=1:total_steps
r_funct = r_funcs(r_funct, t, 0, 0, 0, 0, 4);
end

plot(r_funct(2,:),r_funct(1,:))

%creates area plot
area_matrix = [];
for t=1:length(r_funct)
	area = causal_nonlocal_int(t,r_funct,e,step,rect_thickness);
	area_matrix = [area_matrix area];
end

plot(r_funct(2,:),area_matrix,'Color','g')
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('Area $(\frac{1}{t^2}$','FontSize',14,'interpreter','latex');
title('$\epsilon$ Dependence','FontSize',18,'FontWeight','bold','interpreter','latex');
hold on;

axis([-10 0 -30 1])

%creates analytic comparison plot
r_funcl = deal(NaN(2,total_steps));

for i=1:total_steps
    r_funcl(2,i) = deal(t0 + (i-1)*step);
end

%t
%r_funcl(1,:) = r_funcl(2,:).*(log(-1*(r_funcl(2,1)-r_funcl(2,:)))-1) + t0; 

%1/t^2 old
%r_funcl(1,:) = 1./(r_funcl(2,:)*r_funcl(2,1)) - 1./r_funcl(2,:).^2 + (1./r_funcl(2,:).^2).*log(abs((r_funcl(2,:)-r_funcl(2,1))/r_funcl(2,1)));

%1/t^2 new
r_funcl(1,:) = 1./(r_funcl(2,:)*r_funcl(2,1)) - 1./r_funcl(2,:).^2 + 1./r_funcl(2,:).^2.*log(Mr*r_funcl(2,:).*(r_funcl(2,:)-r_funcl(2,1))/r_funcl(2,1));
plot(r_funcl(2,:),r_funcl(1,:),'Color', 'r')
legend({'.01','.001','.0005','Analytic'},'Position',[.1,.1,.3,.2])

diffs = r_funcl(1,:) - area_matrix;

plot(r_funct(2,:),diffs,'Color','b')
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('Analytic - Numerical','FontSize',14,'interpreter','latex');
title('Residuals for $t$','FontSize',18,'FontWeight','bold','interpreter','latex');


%NUMERICAL ERROR PLOTS

%analytical derivatives

%t array initialization
for i=1:total_steps
    t(i) = deal(t0 + (i-1)*step);
end

a = (-t).^(2/3)/10^(2/3);
a_dot = -2/3/10^(2/3)*(-t).^(-1/3);
a_double_dot = -2/9/10^(2/3)*(-t).^(-4/3);
a_triple_dot = -8/27/10^(2/3)*(-t).^(-7/3);

a_residuals = scale_factor(1,4:end)-a;
a_dot_residuals = scale_1deriv(1,4:end)-a_dot;
a_double_dot_residuals = scale_2deriv(1,4:end)-a_double_dot;
a_triple_dot_residuals = scale_3deriv(1,4:end)-a_triple_dot;

plot(t,a_triple_dot_residuals);
xlabel('Time (s)','FontSize',14,'interpreter','latex');
ylabel('Numerical-Analytic','FontSize',14,'interpreter','latex');
title('$a(t)$ Third Deriv Residuals','FontSize',18,'FontWeight','bold','interpreter','latex');

