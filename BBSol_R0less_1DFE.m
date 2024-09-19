clear all
clc


params = [1.25382975108040	8.08568819620768	0.578027454719331	6.83970198884737e-07];
params = [1.25382975108040	1.908568819620768	0.578027454719331	6.83970198884737e-07 0.09234664];%gives R0<1 and stability DFE. when second param \bets_d is reduced to 1.98568819620768 and 
% down from estimated 8.08568819620768


[ I,Iv,S,Sv,R,B,Bd,R0, time_mesh] = Epi_Zika(params);

R0
figure(1)
plot(time_mesh,I, '-b','LineWidth',3)
% legend('Infected  when R0 <1 ','Infected when R0>1')
xlabel('Days','FontSize',14,'FontName','Sans-serif' );
ylabel('Infected indiviual','FontSize',14,'FontName','Sans-serif' );
ax = gca;
ax.YAxis.Exponent = 0;
set(gca,'LineWidth',3,'FontSize',14,'FontName','Sans-serif');
 title('(b)')
 ax = gca;
ax.TitleHorizontalAlignment = 'left';

%hold off





function [ I,Iv,S,Sv,R,B,Bd,R0, time_mesh] = Epi_Zika(params)




T = 300;%500 % Final Time
Tau = 10; %Final Infectiousness Time
dt = 0.1; %time step=dt

time_mesh = 0:dt:T;
tau_mesh = 0:dt:Tau;




mu_h = 1/(79*365);% per day % host natural death rate
Lambda_h = 20612439*mu_h;%19380000*mu_h; per day % host recruitment rate 


%betav =  0.093007119888124; % transmission rate from infected mosquito to susceptible host class


N = length(time_mesh);
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


% keep your reproduction number formula here

Sv = zeros(N,1);
Iv = zeros(N,1);%00

S = zeros(N,1);

% infected = zeros(N,M_NP);
inf_old = zeros(1,M);
inf_new = zeros(1,M);

I = zeros(N,1);
%C = zeros(N,1);
R = zeros(N,1);
NP = zeros(N,1);

%C(1) = 0;
Nv =  55000;%500;
Lambda_v = Nv*1/10; % 0.98*Nv*1/10;
 muv = 1/10;%1/30;  muv = 1/10;%1/30;  
Sv(1) = 0.98*Nv;%50000;
Iv(1) = 0.02*Nv;%3000;

S(1) =  1000;
inf_old(:) =0.3;
I(1) = sum(inf_old(1:end-1).*dt); 
R(1) = 0;
NP(1) = S(1) + I(1) + R(1);


 for n = 1:N-1
     
     id_cases = (sum(beta_d(1:end-1).*inf_old(1:end-1).*dt))/NP(n);
    
     i_cases = (sum(beta(1:end-1).*inf_old(1:end-1).*dt))/NP(n);
     
     
      
     Sv(n+1) = (Lambda_v*dt + Sv(n))/(1 + i_cases*dt + muv*dt);
           
     Iv(n+1) = (Iv(n) + Sv(n+1)*i_cases*dt)/(1 + muv*dt); 
    
     S(n+1) = (Lambda_h*dt + S(n))/( 1 + id_cases*dt + betav*Iv(n+1)*dt/NP(n) + mu_h*dt);
    
     for k=1:M-1
      
        inf_new(k+1) = (inf_old(k))/(1 + dt*(gamma(k) + alpha(k) + mu_h));

     end
    

    inf_new(1) = betav*S(n+1)*Iv(n+1)/NP(n) + S(n)*id_cases;%boundary condtion
    
   %  C(n+1) = C(n) + inf_new(1);

     I(n+1) = sum(inf_new(1:end-1).*dt);
 
     
     R(n+1) = (R(n) + sum(gamma(1:end-1).*inf_new(1:end-1).*dt))/(1 + mu_h*dt);
     
     NP(n+1) = S(n+1) + I(n+1) + R(n+1); 
    
     inf_old = inf_new;

     
 end

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
end
