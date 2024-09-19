
function bif4

clear all
close all
clc
format long


B=.9;
Lambda = .6; 
%mu =.1;

mu =.3;
Lambdav=2;
betav = 3;
D = 0.01;
%   betav = 5;
 

muv=1.5;
Bd = @(R0) (R0 - (betav*Lambdav*mu*B)/(Lambda*muv^2));
a0 = @(R0) Lambda*muv*(1-D)*(muv*(1-D) + B*mu);
b0 = @(R0) Lambda*muv*mu*(muv*(1 - R0) + muv*(1-D) + mu*B*(1-Bd(R0)) ...
      - muv*(1-Bd(R0)));
c0 = @(R0) Lambda*((muv*mu)^2)*(1-R0);

fcontour(@(R0,X) a0(R0)*X^2 + b0(R0)*X + c0(R0),[0 2 0 1],'-k','LineWidth',3,'LevelList',0)


hold on
 xlabel('$\mathcal{R}_{0}$','Interpreter','latex','FontSize',20,'FontName','Sans-serif')
 set(gca,'LineWidth',2,'FontSize',8,'FontName','Sans-serif');
 ylabel('X','FontSize',14,'FontName','Sans-serif')
  title({'Backward Bifurcation with Endemic Equilibria', 'when $\mathcal{R}_{0} < 1$'},'Interpreter','latex','FontSize',15,'FontName','Sans-serif','FontWeight','normal')

muv=.85;
Bd = @(R0) (R0 - (betav*Lambdav*mu*B)/(Lambda*muv^2));
a0 = @(R0) Lambda*muv*(1-D)*(muv*(1-D) + B*mu);
b0 = @(R0) Lambda*muv*mu*(muv*(1 - R0) + muv*(1-D) + mu*B*(1-Bd(R0))...
     - muv*(1-Bd(R0)));
c0 = @(R0) Lambda*((muv*mu)^2)*(1-R0);
fcontour(@(R0,X) a0(R0)*X^2 + b0(R0)*X + c0(R0),[0 2 0 1],'-r','LineWidth',3,'LevelList',0)
muv=.8;
Bd = @(R0) (R0 - (betav*Lambdav*mu*B)/(Lambda*muv^2));
a0 = @(R0) Lambda*muv*(1-D)*(muv*(1-D) + B*mu);
b0 = @(R0) Lambda*muv*mu*(muv*(1 - R0) + muv*(1-D) + mu*B*(1-Bd(R0)) ...
     - muv*(1-Bd(R0)));
c0 = @(R0) Lambda*((muv*mu)^2)*(1-R0);
fcontour(@(R0,X) a0(R0)*X^2 + b0(R0)*X + c0(R0),[0 2 0 1],'-b','LineWidth',3,'LevelList',0)
% 
muv=.75;
Bd = @(R0) (R0 - (betav*Lambdav*mu*B)/(Lambda*muv^2));
a0 = @(R0) Lambda*muv*(1-D)*(muv*(1-D) + B*mu);
b0 = @(R0) Lambda*muv*mu*(muv*(1 - R0) + muv*(1-D) + mu*B*(1-Bd(R0)) ...
      - muv*(1-Bd(R0)));
c0 = @(R0) Lambda*((muv*mu)^2)*(1-R0);
fcontour(@(R0,X) a0(R0)*X^2 + b0(R0)*X + c0(R0),[0 2 0 1],'-g','LineWidth',3,'LevelList',0)
% 
muv=.7;
Bd = @(R0) (R0 - (betav*Lambdav*mu*B)/(Lambda*muv^2));
a0 = @(R0) Lambda*muv*(1-D)*(muv*(1-D) + B*mu);
b0 = @(R0) Lambda*muv*mu*(muv*(1 - R0) + muv*(1-D) + mu*B*(1-Bd(R0))...
      - muv*(1-Bd(R0)));
c0 = @(R0) Lambda*((muv*mu)^2)*(1-R0);
fcontour(@(R0,X) a0(R0)*X^2 + b0(R0)*X + c0(R0),[0 2 0 1],'-m','LineWidth',3,'LevelList',0)
legend('muv=1.5','muv=0.85','muv=0.8','muv=0.75','muv=.7')
% 
% 


end

