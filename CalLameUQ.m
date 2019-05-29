clear; set(0,'defaultaxesfontsize',20); format long
%%% p2.m smoothing problem for the deterministic logistic map (Ex. 1.4)
%% setup

%J=20;% number of steps
gamma=1e-2;% observational noise variance is gamma^2
C0=1e12;% prior initial condition variance
m0=5e6;% prior initial condition mean
sd=10;rng(sd);% choose random number seed

Info.rhoS = 1.1;
Info.mu = 5.75e6;
Info.lambda = 1.7e6;% dynamics determined by Info

dt = 1e-2;
N = 1;
J = 20;

%% truth
data_t = CalLameParameter(dt, N, J, Info);
[~, dim] = size(data_t);
y = zeros(size(data_t));
for j = 1:J
    y(j, :) = data_t(j, :) + gamma * randn(1, dim);
end

%% solution

mu0 = 0.1e6:0.005e6:9.9e6;% construct vector of different initial data
Phidet=zeros(length(mu0),1);Idet=Phidet;Jdet=Phidet;% preallocate space
vv=zeros(J, dim);% preallocate space
% loop through initial conditions vv0, and compute log posterior I0(vv0)
for j=1:length(mu0)
    fprintf('process: %d/%d\n', j, length(mu0))
    Info.mu = mu0(j);
    Jdet(j)=1/2/C0 * (mu0(j)-m0)^2;% background penalization
    data = CalLameParameter(dt, N, J, Info);
    for i=1:J
        Phidet(j)=Phidet(j)+1/2/gamma^2 * norm(data(i,:) - y(i,:), 2)^2;% misfit functional
    end
    Idet(j)=Phidet(j)+Jdet(j);
end

constant=trapz(mu0,exp(-Idet));% approximate normalizing constant
P=exp(-Idet)/constant;% normalize posterior distribution
prior=normpdf(mu0,m0,sqrt(C0));% calculate prior distribution

figure(2),hold on,
plot(mu0,prior,'r','LineWidth',2)
plot(mu0,P,'r--','LineWidth',2)
plot([5.75e6,5.75e6],[0,max(P)],'b-.','LineWidth',2)
xlabel '\mu_0',
legend 'prior' J=20 'Truth'
hold off