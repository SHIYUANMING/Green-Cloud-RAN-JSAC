function [H]=channel_realization(Area,L,K,N_set)
%%Generate Channel Model

%INPUT:
%L      =# RRHs
%K      =# MUs
%N_set    =# antennas at each RRH
%Area      =length of the area

%OUTPUT:
% H     = (sum(N_set) x K) channel matrix
%%%%%%%%%%%%%%%%%%%%%Network Size%%%%%%%%%%%%%%%%%%%%
M_position=Area*(rand(2,K)-0.5);  %%Position of MUs
R_position=Area*(rand(2,L)-0.5);  %%Position of RAUs

%%%%%%%%%%%%%%%%%%%%Large Scale Fading%%%%%%%%%%%%%%%
for l=1:L
    for k=1:K
       d=norm(R_position(:,l)-M_position(:,k))+10;
       D(l,k)=5.6234*10^(5.355)/(d^(1.88))*sqrt(exp(normrnd(0,6.3)));
    end
end

%%%%%%%%%%%%%%%%%%%%Small Scale Fading%%%%%%%%%%%%%%%

for k=1:1:K
    for l=1:1:L
        if l==1
              H(1:N_set(1),k)=D(l,k)*(normrnd(0,1/sqrt(2),N_set(1),1)+1i*normrnd(0,1/sqrt(2),N_set(1),1));  %%%nosie normalized to 1
        else
         H(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),k)=D(l,k)*(normrnd(0,1/sqrt(2),N_set(l),1)+1i*normrnd(0,1/sqrt(2),N_set(l),1)); 
        end%%%nosie normalized to 1
     end
end   %%Normalized Distribution zero mean and unit variance
        

end


