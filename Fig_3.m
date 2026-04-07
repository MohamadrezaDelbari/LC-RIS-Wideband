%% Impact of N and frequency and receive power
clear variables
clc
close all
% Delete all the path cashes
restoredefaultpath
addpath('./functions')



%% Parameters
Powerdbm=10;                     % Power in dB
Power=db2pow(Powerdbm-30);       % Power
p_bs=[10 10 5];
p_irs=[0 0 0];
p_mu=[0 0 0];
Krice_dB=[10,10,10];
[K,~] = size(p_mu);
SNR_thr(1:K)=10;
h_blk=db2pow(-40);
f=60*10^9;                       % Frequency
c=3*10^8;                        % Speed of the light in vacuum
lambda=c/f;                      % Wave length
kappa=2*pi/lambda;               % Wave number
beta=((lambda/(4*pi))^2);
d_bs_x=lambda/2;                 % Element distance in BS in x-axis
d_bs_z=lambda/2;                 % Element distance in BS in z-axis
d_irs_y=lambda/2;                % Element distance in RIS in y-axis
d_irs_z=lambda/2;                % Element distance in RIS in z-axis
d0=1;
N0=db2pow(-174)/1000;
W_total=8.64*10^9;
W=W_total/2048;
Nf=db2pow(6);
var_noise_base=N0*W*Nf;
norm_factor = var_noise_base;
var_noise = var_noise_base/norm_factor;
factor_Hi = 1/sqrt(sqrt(norm_factor));
factor_Hr = factor_Hi;
factor_Hd = 1/sqrt(norm_factor);
eta=[2 2 2];                     % Path-loss exponent
N_mu=1;                          % Number of MU antennas
N_y=100;                         % Number of RIS elements in y-axis
N_z=1;                           % Number of RIS elements in z-axis
S=10;                            % Number of scatterers
% Scatterers
Vscatter{1}=[p_bs;mean(p_mu,1)]; % direct channel
Vscatter{2}=[p_irs;p_bs]; % BS-RIS
Vscatter{3}=[p_irs;mean(p_mu,1)]; % RIS-MU

Nend=100;
for N=2:Nend
   N_bs_x=N;
   N_bs_z=1;
   i=1;
   % Design precoder for center frequency
    ParamC = struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
           'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
           'eta',eta,'d0',d0,'beta',beta,'h_blk',h_blk,...
           'lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'N_mu',N_mu,...
           'N_y',N_y,'N_z',N_z,'kappa',kappa,...
           'd_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
           'Krice_dB',Krice_dB,'Vscatter',Vscatter,'S',S);
       [H_d,H_i,H_r,Param_output] = func_channel(ParamC);
       a_bs_mu_NF=Param_output.a_bs_mu_NF;
       H_d_los_NF=Param_output.H_d_los_NF;
       SNR_base=abs(H_d_los_NF*a_bs_mu_NF(:,1)/N)^2*Power;
   for kappai=2*pi/c*(56e9:0.1e9:64e9)
   %for kappai=2*pi/c*(60e9:1e9:64e9)
   ParamC = struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
           'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
           'eta',eta,'d0',d0,'beta',beta,'h_blk',h_blk,...
           'lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'N_mu',N_mu,...
           'N_y',N_y,'N_z',N_z,'kappa',kappai,...
           'd_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
           'Krice_dB',Krice_dB,'Vscatter',Vscatter,'S',S);
       [H_d,H_i,H_r,Param_output] = func_channel(ParamC);
       H_d_los_NF=Param_output.H_d_los_NF;
       SNR(N-1,i)=pow2db(abs(H_d_los_NF*a_bs_mu_NF(:,1)/N)^2*Power/SNR_base);
       i=i+1;
   end
end
minSNR=min(SNR,[],2);
N=2:Nend;
figure
fi=56e9:0.1e9:64e9;
plot(fi,SNR(7,:),fi,SNR(19,:),fi,SNR(59,:),fi,SNR(99,:))
figure
plot(N,minSNR)
