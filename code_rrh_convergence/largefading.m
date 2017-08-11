function [D]=largefading(Area,L,K,N_set)
%%Generate Channel Model

%INPUT:
%L      =# RRHs
%K      =# MUs
%N_set    =# antennas at each RRH
%Area      =length of the area

%OUTPUT:
% H     = (sum(N_set) x K) channel matrix
%%%%%%%%%%%%%%%%%%%%%Network Size%%%%%%%%%%%%%%%%%%%%
% M_position=Area*(rand(2,K)-0.5);  %%Position of MUs
% R_position=Area*(rand(2,L)-0.5);  %%Position of RAUs

%% Meshgrid position
%M_position=meshgrid( linspace(-Area/2,Area/2,K), linspace(-Area/2,Area/2,K));  %%Position of MUs
% M_position=Area*(rand(2,K)-0.5);  %%Position of MUs
% 
% R_position=zeros(2,L);
% [X_temp, Y_temp]=meshgrid( linspace(-Area/2,Area/2,sqrt(L)), linspace(-Area/2,Area/2,sqrt(L)));  %%Position of RAUs
% X_vec=vec(X_temp); Y_vec=vec(Y_temp);
% for l=1:L
%     R_position(1,l)=X_vec(l);
%     R_position(2,l)=Y_vec(l);
% end

%%%%%%%%%%%%%%%%%%%%Large Scale Fading%%%%%%%%%%%%%%%
for l=1:L
    for k=1:K
        %d=norm(R_position(:,l)-M_position(:,k))+10;
        %D(l,k)=5.6234*10^(5.355)/(d^(1.88))*sqrt(exp(normrnd(0,6.3)));%%Noise power -102dBm
        %D(l,k)=5.6234*10^(5.355)/(d^(1.88));%%Noise power -102dBm
        D(l,k)=1;
    end
end

for k=1:K
    temp=randperm(L);
    
%     D(temp([1:L/2]),k)=ones(L/2,1);
%     D(temp([L/2+1:L]),k)=0.5*ones(L/2,1);

%     D(temp([1:L/3]),k)=ones(L/3,1);
%     D(temp([L/3+1:2*(L/3)]),k)=0.7*ones(L/3,1);
%     D([2*(L/3)+1:L],k)=0.5*ones(L/3,1);
end
    
end