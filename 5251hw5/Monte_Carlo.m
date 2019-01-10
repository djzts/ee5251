clc
E_t=[];Error_t=[];mu_t=[];
qlist=[0.1 0.5 1 10 20 100 1000 1e4 1e5];
count=100;

%% body
for z=1:length(qlist)
q_nom=qlist(z);
    for i=1:count
        example4_2M;
        E_t=[E_t;e_t];
        Error_t=[Error_t,error500];
        mu_t=[mu_t;mu];
    end
    [m1,n1]=size(Error_t);
    rho_t=[];
    for k=1:count
        for j=1:count
            part1=(sum(Error_t(k,:)*Error_t(k,:)')*sum(Error_t(j,:)* Error_t(j,:)'))^(-0.5);
            part2=sum(Error_t(k,:)*Error_t(j,:)');
            rho_t=[rho_t,1/(sqrt(count))*part1*part2];
        end
    end
    q_nom
    epsilon=mean(E_t);
    rho_t=reshape(rho_t,[count,count]);
    mu_k=mean(mu_t);
    muk1=sum(mu_k<=1.96/10 & mu_k>=-1.96/10)/500;
    rhot1=sum(sum(rho_t<=1.96/10 & rho_t>=-1.96/10))/(count^2);

%% Result
    disp([mean(mu_k),muk1,rhot1]);
%disp(rho_t)
end

disp("done");
