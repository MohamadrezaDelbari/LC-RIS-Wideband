%% Assume we have the matrix for different frequencies

%%{
clear variables
clc
close all
restoredefaultpath
addpath('./functions')
addpath(genpath('/home/mdbab3f/Documents/MATLAB/Me/cvx'), '-begin');


load('Data/VV.mat');
v=VV(:,2:5);
K=4;
N=30;
%v=rand(N,K)*10;
%v(3,:)=0;
normalizedpower=db2pow(2*[-1.5 0 -0.4 -2]); % Main 2 is for amplitude to power
%normalizedpower=normalizedpower.^2;
%normalizedpower=db2pow([-1.5 0 -0.4 -2]); % Fake

for k=1:K
    a1(:,k)=exp(1i*voltagephase(v(:,k))/180*pi)';
    a1(:,k)=a1(:,k)/a1(3,k);
    A1(:,:,k)=a1(:,k)*a1(:,k)';

    f(k)=60e9+(k-floor(K/2))*2e9;
    beta(k)=1+2.4*(f(k)/60e9-1);
    a2(:,k)=a1(:,k).^beta(k);
    A2(:,:,k)=a2(:,k)*a2(:,k)';
end

%%}

%% Only considering the channel impact not element

A=A1;

%S_0=A(:,:,floor(K/2))';
S_0=A(:,:,K)';
Imax=12;
eta=0.01;


for ii=1:Imax

[V_2,D_2]=eig(S_0);
[~,index_v]=max(diag(D_2));
V_2=V_2(:,index_v);

cvx_begin
%cvx_solver MOSEK
variable S(N,N) complex hermitian; 
variable gamma1
maximize gamma1-eta*(real(trace(S)-trace(V_2*V_2'*(S-S_0))))
subject to
diag(S)==1; 
S== hermitian_semidefinite(N);

for k=1:K
    gamma1<=real(trace(A(:,:,k)*S))*normalizedpower(k);
end
cvx_end

eta=min(5*eta,5*10^5);
S_0=S;

end

s=(sqrt(max(max(D_2)))*V_2)';
s=s/s(3);
vs=phasevoltage(mod(angle(s),2*pi)*180/pi).';


%% With considering element frequency response
A=A2;
S_0=S;
Imax=12;
eta=0.01;


for ii=1:Imax

[V_2,D_2]=eig(S_0);
[~,index_v]=max(diag(D_2));
V_2=V_2(:,index_v);

cvx_begin
%cvx_solver MOSEK
variable S(N,N) complex hermitian; 
variable gamma2
maximize gamma2-eta*(real(trace(S)-trace(V_2*V_2'*(S-S_0))))
subject to
diag(S)==1; 
S== hermitian_semidefinite(N);

for k=1:K
    gamma2<=real(trace(A(:,:,k)*(S_0.^beta(k)+beta(k)*S_0.^(beta(k)-1).*(S-S_0))))*normalizedpower(k);
end
cvx_end

eta=min(5*eta,5*10^5);
S_0=S;

end

s2=(sqrt(max(max(D_2)))*V_2)';
s2=s2/s2(3);
vs2=phasevoltage(mod(angle(s2),2*pi)*180/pi).';