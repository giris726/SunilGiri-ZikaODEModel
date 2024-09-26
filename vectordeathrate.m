B=.9;
Lambda = .6;

mu =.3;
Lambdav=2;
betav = 3;
D = 0.1;

muv=0.28;% for R0 in BB, muv value less or equal than 0.28 the dicriminant is negative and hence no Endemic
R0=0.9705;

Bd = (R0 - (betav*Lambdav*mu*B)/(Lambda*muv^2));
a0 = Lambda*muv*(1-D)*(muv*(1-D) + B*mu);
b0 = Lambda*muv*mu*(muv*(1 - R0) + muv*(1-D) + mu*B*(1-Bd) ...
      - muv*(1-Bd));
c0 = Lambda*((muv*mu)^2)*(1-R0);

b0*b0-4*a0*c0 %for R0 in BB, muv value less or equal than 0.28 the dicriminant is negative and hence no Endemic