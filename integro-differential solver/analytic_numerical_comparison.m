
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
