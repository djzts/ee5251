clc

% Initialize
dt=0.01;tf=10;t=[0:dt:tf]';m=length(t);
x0=[1;1]; e_t=[0]; Pk=[]; mu=[0];

% System Matricess
a=[0.9999 0.0099;-0.0296 0.9703];b=[0;0];h=[1 0];
r=0.01;q=10;
%phi=c2d(a,b,dt);
phi=a;
upsilon=[0;0.01];

% Uncomment These Lines to Generate Your Own Measurements
x=zeros(m,2);x(1,:)=x0';
y=zeros(m,1);y(1)=h*x(1,:)';
for i = 1:m-1
 x(i+1,:)=(phi*x(i,:)'+upsilon*sqrt(q)*randn(1))';
 y(i+1)=h*x(i,:)';
end
ym=y+sqrt(r)*randn(m,1);




% Trial Values
q_trial=[0.1;0.5;1;10;20;100;1000;1e4;1e5];

% Main Loop 
for j_trial=1:length(q_trial)

% Calculate Gain
q_nom=q_trial(j_trial);
p=dare(phi',h',upsilon*q_nom*upsilon',r);
gain=p*h'*inv(h*p*h'+r);

% Residual
xe=zeros(m,2);xe(1,:)=x0';
cov_res=inv(h*p*h'+r);
e_time=0;j=0;

for i = 1:m-1
 xe(i+1,:)=(phi*xe(i,:)'+phi*gain*(ym(i)-h*xe(i,:)'))';
    if i > 500
        j = j+1;
        e_t(j)=(ym(i+1)-xe(i+1,1))^2*cov_res;
        mu(j)=(ym(i+1)-xe(i+1,1))/sqrt(cov_res);
        e_time=e_time+(ym(i+1)-xe(i+1,1))^2*cov_res;
       
    end
end

e_time=e_time/j;

% Whiteness Test
error500=ym(501:1001)-xe(501:1001,1);
[rho,limit95]=whiteness_test(error500,2);
q_nom
muk=sum(mu<=1.96/sqrt(500) & mu>=-1.96/sqrt(500))/500;

%% result
disp([mean(mu),rho,limit95,muk]);

end
%disp('done');


