% play with parameter values mainly with transmission to get different sensititvity inde
Lambda_v = 1/10;%1/10;%10000; % 
 muv = 1/10;
 mu_h = 1/(79*365); % host natural death rate
Lambda_h = 20612439*mu_h;

params = [1.25382975108040	1.908568819620768	0.578027454719331	6.83970198884737e-07 0.09234664];


dt = 0.1; %time step=dt

Tau = 10;
tau_mesh = 0:dt:Tau;
M = length(tau_mesh);

V = normpdf(tau_mesh,5,3);%Gaussian

%
 beta=zeros(1,M);
beta(:) = params(1).*V'; %taking gaussian value 
 beta_d = zeros(1,M);
beta_d(:) = params(2).*V'; 
alpha = zeros(1,M);
alpha(:) = params(3);%taking all  equal value
gamma =zeros(1,M);
gamma(:) = params(4);%taking all  equal value
betav = params(5);



 s_mesh= 0:0.001:1;
 J= length(s_mesh);
    mu= zeros(1,J);
aalpha =zeros(1,J);
ggamma = zeros(1,J);
mu(:) = mu_h;
 aalpha(:)= params(3);
 ggamma(:) = params(4);
 tsum= ggamma(:) + aalpha(:) + mu(:);
pi= zeros(1,M);
pi(:) = exp(- sum(tsum(1: end-1).*0.001));
%pi(:) = exp(- (params(3)+params(4) +mu_h));
     Bd   = sum(beta_d(1:end-1).*pi(1:end-1).*dt);
     B     = sum(beta(1:end-1).*pi(1:end-1).*dt);

     %Reproduction number
     R0 = (((betav*Lambda_v*mu_h)/((muv^2)*Lambda_h))*B) + Bd;

  S_beta =  (1.25382975108040/R0)*((Lambda_v*betav*mu_h)/(Lambda_h*muv^2))*sum(pi(1:end-1).*dt);

% S_B = (B*Lambda_v*betav*mu_h)/(Lambda_h*muv^2*(Bd + (B*Lambda_v*betav*mu_h)/(Lambda_h*muv^2)));

S_betad = (1.908568819620768/R0)*sum(pi(1:end-1).*dt);

% S_Bd = Bd/(Bd + (B*Lambda_v*betav*mu_h)/(Lambda_h*muv^2)) ;

 S_betav = (betav/R0)*((Lambda_v*mu_h)/(Lambda_h*muv^2))*(1.25382975108040)*sum(pi(1:end-1).*dt);

 S_muv = (-2/R0)*((betav*Lambda_v*mu_h)/(Lambda_h*muv^2))*(1.25382975108040)*sum(pi(1:end-1).*dt);

 S_Lambda_v = (Lambda_v/R0)*((betav*mu_h)/(Lambda_h*muv^2))*(1.25382975108040)*sum(pi(1:end-1).*dt);

 S_Lambda_h = (-1/R0)*((betav*Lambda_v*mu_h)/(Lambda_h*muv^2))*(1.25382975108040)*sum(pi(1:end-1).*dt);

 S_mu_h = (mu_h/R0)*((betav*Lambda_v)/(Lambda_h*muv^2))*(1.25382975108040)*sum(pi(1:end-1).*dt);


 S_value_params = [S_beta, S_betad, S_betav, S_muv, S_Lambda_v, S_Lambda_h, S_mu_h];

 c = categorical({'\beta(\tau)', '\beta_d(\tau)', '\beta_{v}', '\mu_{v}', '\Lambda_{v}', '\Lambda', '\mu'});



figure
bar(c,S_value_params,1)

 

ylabel('Sensitivity Index, S_{\alpha}','FontSize',4)
xlabel('Parameters','FontSize',14)
%ylim([-0.000005 0.000005]) 
ylim([-0.000001 0.000005]) 

set(gca,'LineWidth',3,'FontSize',8)



 

 

 




