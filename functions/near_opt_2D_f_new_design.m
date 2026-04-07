function [S,s] = near_opt_2D_f_new_design(Paramphase,ParamC)
lambda=Paramphase.lambda;
p_bs=Paramphase.p_bs;
pris=Paramphase.pris;
pp_mu=Paramphase.pp_mu;
pp_e=Paramphase.pp_e;
k=2*pi/lambda;
step=Paramphase.step;
N=length(pris);
S_0=Paramphase.S_0;
Power=Paramphase.Power;
kappa=ParamC.kappa;
kappa2=Paramphase.kappa;
gamma=Paramphase.gamma;





Imax=8;

eta=0.00017;

%______________________Near field________________________

prx=pp_mu(:,1);
pry=pp_mu(:,2);
ppr_x=linspace(prx(1),prx(2),(prx(2)-prx(1))/step+1);
ppr_y=linspace(pry(1),pry(2),(pry(2)-pry(1))/step+1);

prxe=pp_e(:,1);
prye=pp_e(:,2);
ppr_xe=linspace(prxe(1),prxe(2),(prxe(2)-prxe(1))/step+1);
ppr_ye=linspace(prye(1),prye(2),(prye(2)-prye(1))/step+1);


for yy=1:length(ppr_y)
    for xx=1:length(ppr_x)
        p_test = [ppr_x(xx) ppr_y(yy) pp_mu(1,3)];
        for i=1:length(ParamC)
            ParamC(i).p_mu=p_test;
        end
        for t=1:length(kappa2)
            kappa_temp=kappa2(t);
            factor(t)=1+2.4*(kappa_temp-kappa)/kappa;
            %factor(t)=kappa_temp/kappa;
            for i=1:length(ParamC)
                ParamC(i).kappa=kappa_temp;
            end
            [H_d,H_i,H_r,Param_output] = func_channel(ParamC);
            v_bs=Param_output(1).a_bs_irs;
            a=diag(H_r{1})*H_i*v_bs(:,1)*sqrt(Power);
            A(:,:,xx,yy,t)=a*a';
        end
    end
end

for yy=1:length(ppr_ye)
    for xx=1:length(ppr_xe)
        p_test = [ppr_xe(xx) ppr_ye(yy) pp_e(1,3)];
        for i=1:length(ParamC)
            ParamC(i).p_mu=p_test;
        end
        for t=1:length(kappa2)
            kappa_temp=kappa2(t);
            for i=1:length(ParamC)
                ParamC(i).kappa=kappa_temp;
            end
            [H_d,H_i,H_r,Param_output] = func_channel(ParamC);
            v_bs=Param_output(1).a_bs_irs;
            ae2=diag(H_r{1})*H_i*v_bs(:,1)*sqrt(Power);
            Ae(:,:,xx,yy,t)=ae2*ae2';
        end
    end
end


%for i_S=1:10
%eta=0.0018;
%vv=rand(1,N)*2*pi;
%v_int=exp(1i*vv);
%S_0=v_int'*v_int;

Rankcost(1:Imax,1:2)=zeros(Imax,2);

for iter=1:2
 t=1;
for yye=1:length(ppr_ye)
for xxe=1:length(ppr_xe)
for yy=1:length(ppr_y)
for xx=1:length(ppr_x)
    for tt=1:length(kappa2)
        A_final(:,t,tt)=reshape((A(:,:,xx,yy,tt)-gamma*Ae(:,:,xxe,yye,tt)).',[],1);
    end
t=t+1;
end
end
end
end
T=t-1;

for ii=1:Imax

[V_2,D_2]=eig(S_0);
[~,index_v]=max(diag(D_2));
V_2=V_2(:,index_v);

cvx_begin
%cvx_solver MOSEK
variable S(N,N) complex hermitian; 
variable gamma
maximize gamma-eta*(real(trace(S)-trace(V_2*V_2'*(S-S_0))))
subject to
diag(S)==1; 
S== hermitian_semidefinite(N);

for tt=1:length(kappa2)
    0<=1-gamma+real(A_final(:,:,tt).'*reshape(S_0.^factor(tt)+factor(tt)*S_0.^(factor(tt)-1).*(S-S_0),[],1))
end
cvx_end
% if norm(S-S_0)<0.01
%     break
% end
eta=min(5*eta,5*10^5);
Rankcost(ii,iter)=real(trace(S_0)-trace(V_2*V_2'*(S_0)));
S_0=S;

end

end
%     Rankfinal(i_S,:)=Rankcost;
% end

%ii=1:Imax;
%AA_rank=[ii',Rankcost(:,1)];
%plot(ii,Rankcost(:,1),ii,Rankcost(:,2))

s=(sqrt(max(max(D_2)))*V_2)';
end