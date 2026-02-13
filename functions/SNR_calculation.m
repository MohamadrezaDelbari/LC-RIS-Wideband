function [SNR,SR] = SNR_calculation(Param,WW,plott)
ii=1;
jj=1;
minu=1000000000;
maxe=-10000000000000;

x_start=Param.x_start;
x_end=Param.x_end;
y_start=Param.y_start;
y_end=Param.y_end;
x_step=Param.x_step;
y_step=Param.y_step;
p_mu=Param.p_mu;
Power=Param.Power;
pp_mu=Param.pp_mu;
pp_e=Param.pp_e;

% Number of points in x and y
Nx = length(x_start:x_step:x_end);
Ny = length(y_start:y_step:y_end);

% Pre-allocate SNR matrix
SNR = zeros(Nx, Ny);  % Adjust dimensions as necessary
ParamC=Param;


% Parallel for-loop over y
for jj = 1:Nx
    xx = x_start + (jj - 1) * x_step;
    for ii = 1:Ny
        yy = y_start + (ii - 1) * y_step;
        p_test = [xx yy p_mu(1,3)];
        for i=1:length(ParamC)
            ParamC(i).p_mu=p_test;
        end

        [H_d,H_i,H_r,Param_output] = func_channel(ParamC);
        v_bs=Param_output(1).a_bs_irs;
        SNR(jj,ii) = pow2db(abs((H_r{1}*WW*H_i + H_d{1,1})*v_bs(:,1))^2*Power);
        if dot(pp_mu(1,:)<=p_test,p_test<=pp_mu(2,:))==3 && db2pow(SNR(jj,ii))<minu
            minu=db2pow(SNR(jj,ii));
        end
        if dot(pp_e(1,:)<=p_test,p_test<=pp_e(2,:))==3 && db2pow(SNR(jj,ii))>maxe
            maxe=db2pow(SNR(jj,ii));
        end
    end
end
SR=log2((1+minu)/(1+maxe));


if plott==1
    x_values = x_start:x_step:x_end;
    y_values = y_start:y_step:y_end;
    [I,J]=meshgrid(y_values,x_values);
    figure;  % Open a new figure for each user
        %SNR_matrix = squeeze(SNR(:,:));  % Extract the matrix for user k
    %plot(y_values,SNR)
    surf(I,J,SNR)
    shading interp
    colorbar('FontSize',15)
    colormap jet
    view(0, 90);
    clim([-20 20])
end

end
