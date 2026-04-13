%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%----- Wideband Illumination Liquid Crystal Reconfigurable Intelligent Surfaces: Modeling, Design, and Experimental Tests
%%----- Please cite: https://arxiv.org/pdf/2604.09214
%%----- Please cite: https://arxiv.org/pdf/2508.04331
%%----------------------------------------------------%%


clear variables
clc
close all
% Delete all the path cashes
restoredefaultpath
addpath('./functions')

%% Parameters
load('Data/Data_100times1f60W8GHzP10')
Powerdbm=10;                     % Power in dB
Power=db2pow(Powerdbm-30);       % Power
f=60*10^9;                       % Frequency
c=3*10^8;                        % Speed of the light in vacuum
lambda=c/f;                      % Wave length
kappa=2*pi/lambda;               % Wave number
N_y=100;                         % Number of RIS elements in y-axis
N_z=1;                           % Number of RIS elements in z-axis
N=N_y*N_z;
N_bs_x=16;                       % Number of BS antennas in x-axis
N_bs_z=16;                       % Number of BS antennas in z-axis
N_mu=1;                          % Number of MU antennas
d_bs_x=lambda/2;                 % Element distance in BS in x-axis
d_bs_z=lambda/2;                 % Element distance in BS in z-axis
d_irs_y=lambda/2;                % Element distance in RIS in y-axis
d_irs_z=lambda/2;                % Element distance in RIS in z-axis
Ly=N_y*lambda/2;                 % y-length in RIS
Lz=N_z*lambda/2;                 % z-length in RIS
D=sqrt(Ly^2+Lz^2);               % RIS diameter
dff=2*D^2/lambda;                % Far-field regime in RIS
eta=[2 2 2];                     % Path-loss exponent
Krice_dB=[-100,10,10];           % K-Ricean factor
S=10;                            % Number of scatterers
fi=56e9:1e9:64e9;
lambdai=c./fi;                   % Wave length
kappai=2*pi./lambdai;            % Wave number
betai=((lambdai/(4*pi)).^2);
f2=56*10^9;
lambda2=c/f2;                    % Wave length
kappa2=2*pi/lambda2;             % Wave number
beta2=((lambda2/(4*pi))^2);
%% Positions
p_bs=[10 10 5];
p_irs=[0 0 0];
p_mu=[6 1 -5];
p_e=[6 -1 -5];
% Scatterers
Vscatter{1}=[p_bs;mean(p_mu,1)]; % direct channel
Vscatter{2}=[p_irs;p_bs]; % BS-RIS
Vscatter{3}=[p_irs;mean(p_mu,1)]; % RIS-MU
%% Normalization factors:
[K,~] = size(p_mu);
SNR_thr(1:K)=10;
h_blk=db2pow(-40);
beta=((lambda/(4*pi))^2);
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

ParamC = struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
           'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
           'eta',eta,'d0',d0,'beta',beta,'h_blk',h_blk,...
           'lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'N_mu',N_mu,...
           'N_y',N_y,'N_z',N_z,'kappa',kappa,...
           'd_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
           'Krice_dB',Krice_dB,'Vscatter',Vscatter,'S',S);
       [H_d,H_i,H_r,Param_output] = func_channel(ParamC);
pris=Param_output.pp_ris_g;

%% Benchmark 1 (All subcarriers, Area)
for kk=1:K
   pp_mu{kk}=[p_mu(kk,:)-1*[1 1 0];p_mu(kk,:)+1*[1 1 0]];
   pp_e=[p_e-1*[1 1 0];p_e+0*[0 1 0]];
   Paramphase = struct('p_bs',p_bs,'pris',pris,'pp_mu',pp_mu,'pp_e',pp_e,...
           'lambda',lambda,'step',1,'Power',Power,'kappa',kappai,'gamma',1000);
   %[Wno1_f,wno1_f]=near_opt_2D_f(Paramphase,ParamC);
   WW_opt_f{kk}=diag(wno1_f);
end

%% Benchmark 2 (Center frequency, Area)
for kk=1:K
   Paramphase = struct('p_bs',p_bs,'pris',pris,'pp_mu',pp_mu,'pp_e',pp_e,...
           'lambda',lambda,'step',1,'Power',Power,'gamma',1000);
   %[Wno1,wno1]=near_opt_2D(Paramphase,ParamC);
   WW_opt{kk}=diag(wno1);
end

%% Benchmark 3 (All subcarriers, Point)
for kk=1:K
   Paramphase = struct('p_bs',p_bs,'pris',pris,'pp_mu',[p_mu;p_mu],'pp_e',[p_e;p_e],...
           'lambda',lambda,'step',1,'Power',Power,'kappa',kappai,'gamma',1000);
   %[Wno1_f_point,wno1_f_point]=near_opt_2D_f(Paramphase,ParamC);
   WW_opt_f_point{kk}=diag(wno1_f_point);
end

%% Proposed (SDP)
for kk=1:K
   Paramphase = struct('p_bs',p_bs,'pris',pris,'pp_mu',pp_mu,'pp_e',pp_e,...
           'lambda',lambda,'step',1,'Power',Power,'kappa',kappai,'gamma',1,'S_0',wno1_f'*wno1_f);
   %[Wno1_f_new_design,wno1_f_new_design]=near_opt_2D_f_new_design(Paramphase,ParamC);
   WW_opt_f_new_design{kk}=diag(wno1_f_new_design);
end

%% Proposed (Scalable)
for kk=1:K
   Paramphase = struct('p_bs',p_bs,'pris',pris,'pp_mu',pp_mu,'pp_e',pp_e,...
           'lambda',lambda,'step',1,'Power',Power,'kappa',kappai,'gamma',0.7,'S_0',wno1_f'*wno1_f,'s_opt',wno1_f_new_design);
   [Wno1_f_scalable,wno1_f_scalable]=near_scalable_2D_f_new_design(Paramphase,ParamC);
   WW_f_scalable{kk}=diag(wno1_f_scalable);
end

%save('Data/Data_100times1f60W8GHzP10','Wno1','wno1','WW_opt_new','Wno1_f','wno1_f','WW_opt_f_new','WW_opt_f_new_design','wno1_f_new_design','wno1_f_point')

%% Calculating the secure rate on a range of frequency %%%%%%%%%%%%%%%%%
clear SR_opt SR_opt_f_new SR_opt_f_new_design
t=1;
for fi=50e9:1e9:70e9
   lambdai=c/fi;                      % Wave length
   kappai=2*pi/lambdai;                   % Wave number
   betai=((lambdai/(4*pi)).^2);
   Param=struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
           'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
           'eta',eta,'d0',d0,'beta',betai,'h_blk',h_blk,...
           'lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'N_mu',N_mu,'kappa',kappai,...
           'N_y',N_y,'N_z',N_z,'d_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
           'Krice_dB',Krice_dB,'Vscatter',Vscatter,'S',S,'x_start',2,'x_end',10, ...
           'y_start',-4,'y_end',4,'x_step',1,'y_step',1,'K',K,'Power',Power,'pp_mu',pp_mu,'pp_e',pp_e);
       factor=1+2.4*(kappai-kappa)/kappa;
       [SNR_opt_f_new(:,:),SR_opt_f_new(t)]=SNR_calculation(Param,WW_opt_f{1}.^factor,0); %Benchmark 1
       [SNR_opt(:,:),SR_opt(t)]=SNR_calculation(Param,WW_opt{1}.^factor,0);  %Benchmark 2
       [SNR_opt_f_point(:,:),SR_opt_f_point(t)]=SNR_calculation(Param,WW_opt_f_point{1}.^factor,0); %Benchmark 3
       [SNR_opt_f_new(:,:),SR_opt_f_new_design(t)]=SNR_calculation(Param,WW_opt_f_new_design{1}.^factor,0); % Proposed
       [SNR_f_scalable(:,:),SR_f_scalable(t)]=SNR_calculation(Param,WW_f_scalable{1}.^factor,0); % Proposed_S
       t=t+1;
end
%% Secrecy rate (bits/symbol) vs frequency (Figure 7)
figure
fi=50e9:1e9:70e9;
fii=50e9:0.1e9:70e9;
SR_opt_int=max(interp1(fi,SR_opt,fii,"spline"),0);
SR_opt_f_new_int=max(interp1(fi,SR_opt_f_new,fii,"spline"),0);
SR_opt_f_point_int=max(interp1(fi,SR_opt_f_point,fii,"spline"),0);
SR_opt_f_new_design_int=max(interp1(fi,SR_opt_f_new_design,fii,"spline"),0);
SR_f_scalable_int=max(interp1(fi,SR_f_scalable,fii,"spline"),0);
plot(fii,SR_opt_int,fii,SR_opt_f_new_int,fii,SR_opt_f_point_int,fii,SR_opt_f_new_design_int,fii,SR_f_scalable_int)
AA=[fii', SR_opt_int',SR_opt_f_new_int',SR_opt_f_new_design_int',SR_opt_f_point_int',SR_f_scalable_int'];

%% Heat maps (Figure 8)
   fi=60e9;                          % Frequency of the signal for heatmap
   lambdai=c/fi;                      % Wave length
   kappai=2*pi/lambdai;                   % Wave number
   betai=((lambdai/(4*pi)).^2);
   Param=struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
           'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
           'eta',eta,'d0',d0,'beta',betai,'h_blk',h_blk,...
           'lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'N_mu',N_mu,'kappa',kappai,...
           'N_y',N_y,'N_z',N_z,'d_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
           'Krice_dB',Krice_dB,'Vscatter',Vscatter,'S',S,'x_start',2,'x_end',10, ...
           'y_start',-4,'y_end',4,'x_step',1,'y_step',1,'K',K,'Power',Power,'pp_mu',pp_mu,'pp_e',pp_e);
       factor=1+2.4*(kappai-kappa)/kappa;
       [SNR_f_scalable(:,:),SR_f_scalable(t)]=SNR_calculation(Param,WW_f_scalable{1}.^factor,1); % Choose which one do you want to plot (Example here: Proposed_S)

hold on

x_rect = [-0.2 2.2 2.2 -0.2]; y_rect = [4.8 4.8 7.2 7.2]; z_rect = ones(1,4) * 30;  % just above the surface
fill3(x_rect, y_rect, z_rect, 'g', ...       % Use 'g' or [0 1 0] for green
   'FaceAlpha', 0, ...                      % Transparent fill
   'EdgeColor', 'g', ...                    % Green edge
   'LineWidth', 5)
x_rect = [-2.2 -0.8 -0.8 -2.2]; y_rect = [4.8 4.8 6.2 6.2]; z_rect = ones(1,4) * 30;  % just above the surface
fill3(x_rect, y_rect, z_rect, 'r', ...      
   'FaceAlpha', 0, ...                      % Transparent fill
   'EdgeColor', 'r', ...                    % Red edge
   'LineWidth', 5)
%title('Benchmark 1','FontSize',28,'FontWeight','bold')
xlabel('$\mathsf{y[m]}$','FontSize',22,'FontWeight','bold','Interpreter','latex')
ylabel('$\mathsf{x[m]}$','FontSize',22,'FontWeight','bold','Interpreter','latex')
exportgraphics(gca, fullfile('Figure', sprintf('Proposed_S_%d.pdf', fi/1e9)), 'ContentType', 'vector')
