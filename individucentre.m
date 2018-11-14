% simulation parameters

%reading emprirical data
data = readtable('MOMv3.3.txt','Delimiter','\t','ReadVariableNames',false,'HeaderLines',0);
empirical = data{:,8};

xmin = 1.8; % lower bound
x0 = 40; % founder body size
n = 5000; % num. species at equilbrium
beta = 1/n; % baseline extinction rate
rho = 0.025; % rate of extinction increase
nu = 1.6; % mean species lifetime (My)

%initializing the parameters for Cope's rule
tau = 60; % total simulation time (My)
c(1) = 0.33; % log-lambda intercept
c(2) = 1.30; % log-size intercept
c(3) = 8; %param (second virage)

delta = 0.04; % systematic bias (Cope’s rule)
sigma = 0.63; % variance
alpha = 0.30; % power-law tail

% data structure set up
tmax = round((tau/nu)*n);
xmax = 10^15;

%initializing the vectors of mass with and without bias
xsb = -Inf*ones(ceil(1.5*n),1); %sans biais
xb = xsb; %avec biais
xbp = xsb; %avec biais "ameliore"
xsb(1) = x0; 
xb(1) = x0;
xbp(1) = x0;

%initializing loop parameters
kdt = 5000;
nsb = 1;
nsSb = 1;
nsbp = 1;
nkb = 0;
nksb = 0;
nkbp = 0;
kd = 1;
f_stop = 0;

% begin main loop
while ~f_stop
% begin cladogenesis step
pairb = [ceil(nsb*rand(1)) nsb+1];
pairsb = [ceil(nsSb*rand(1)) nsSb+1];
pairbp = [ceil(nsbp*rand(1)) nsbp+1]; 
massb = xb(pairb(1),1);
massSb = xsb(pairsb(1),1);
massbp = xbp(pairbp(1),1);

L1b = massb/xmin; % lower bound
L2b = xmax/massb; % upper bound
L1sb = massSb/xmin;
L2sb = xmax/massSb;
L1bp = massbp/xmin;
L2bp = xmax/massbp;

%model of Cope’s rule
%bias plus
if log10(massbp)<c(2)
    % increased bias for small sizes
    mu2 = (-c(1)/c(2))*log10(massbp)+c(1)+delta;
elseif log10(massbp)<c(3)
    % increased bias for small sizes #2
    mu2 = (-c(1)/c(3))*log10(massbp)+c(1)+delta;
else
% uniform bias for large sizes
    mu2 = delta;
end;
%with bias
if log10(massb)<c(2)
    % increased bias for small sizes
    mu = (-c(1)/c(2))*log10(massb)+c(1)+delta;
else
% uniform bias for large sizes
mu = delta;
end;


% Monte Carlo draw of growth factors - with bias
tt1 = [0 0];
while any(tt1<1/L1b | tt1>L2b)
% F(lambda) with power-law tails
tt1 = exp(randn(2,1)*sigma*mu).* ...
((rand(2,1).* ...
(1-1./L1b)+1./L1b).^alpha)./ ...
((rand(2,1).* ...
(1-1./L2b)+1./L2b).^alpha);
end;

% Monte Carlo draw of growth factors - without bias
tt2 = [0 0];
while any(tt2<1/L1sb | tt2>L2sb)
% F(lambda) with power-law tails
tt2 = exp(randn(2,1)*sigma*delta).* ...
((rand(2,1).* ...
(1-1./L1sb)+1./L1sb).^alpha)./ ...
((rand(2,1).* ...
(1-1./L2sb)+1./L2sb).^alpha);
end;

% Monte Carlo draw of growth factors - bias plus
tt3 = [0 0];
while any(tt3<1/L1bp | tt3>L2bp)
% F(lambda) with power-law tails
tt3 = exp(randn(2,1)*sigma*mu2).* ((rand(2,1).* (1-1./L1bp)+1./L1bp).^alpha)./ ...
((rand(2,1).* (1-1./L2bp)+1./L2bp).^alpha);
end;

xb(pairb) = massb.*tt1;
xsb(pairsb) = massSb.*tt2;
xbp(pairbp) = massbp.*tt3;
kd = kd+2;
nsb = nsb+1;
nsSb = nsSb+1;
nsbp = nsbp+1;
% end cladogenesis step
% begin extinction step

% power-law model of extinction risk
%with bias
klb = rand(nsb,1) < ...
10.^(rho*log10(xb(1:nsb))+log10(beta));
if sum(klb)>0
xb(1:sum(~klb)) = xb(~klb);
xb(sum(~klb)+1:nsb) = ...
repmat([-Inf],sum(klb),1);
nsb = sum(~klb);
nkb = nkb+sum(klb);
end;
%without bias
klsb = rand(nsSb,1) < ...
10.^(rho*log10(xsb(1:nsSb))+log10(beta));
if sum(klsb)>0
xsb(1:sum(~klsb)) = xsb(~klsb);
xsb(sum(~klsb)+1:nsSb) = ...
repmat([-Inf],sum(klsb),1);
nsSb = sum(~klsb);
nksb = nksb+sum(klsb);
end;
%bias plus
klbp = rand(nsbp,1) < ...
10.^(rho*log10(xbp(1:nsbp))+log10(beta));
if sum(klbp)>0
xbp(1:sum(~klbp)) = xbp(~klbp);
xbp(sum(~klbp)+1:nsbp) = ...
repmat([-Inf],sum(klbp),1);
nsbp = sum(~klbp);
nkbp = nkbp+sum(klbp);
end;
% end extinction step
% begin check stop-criteria
if kd>=tmax, f_stop = 1; end;
% end check stop-criteria
end;
% end main loop

%get the interesting values in the vectors
xsb = xsb(xsb>xmin & xsb<xmax);
xb = xb(xb>xmin & xb<xmax);
xbp = xbp(xbp>xmin & xbp<xmax);
empirical = empirical(empirical>xmin & empirical < xmax);

%comparing the distributions
mean(log(xsb))
mean(log(xb))
mean(log(xbp))
mean(log(empirical))

%plotting the distributions
[ysb,zsb] = ksdensity(log(xsb));
[yb,zb] = ksdensity(log(xb));
[ybp,zbp] = ksdensity(log(xbp));
[ye,ze] = ksdensity(log(empirical));
figure(1)
plot(zsb,ysb,zb,yb,zbp,ybp,ze,ye);
legend('Sans biais','Avec biais','Biais plus','Empirique')
figure(2)
plot(ze,ye)
