function [w_new] = phase_v_f(w,f,f_r)
v_r=[0 2 4 6 10 20];
w_r=[0 5 50 75 90 100]/100*2*pi;
v_new=interp1(w_r,v_r,w,'linear');
w_new=w*(1+2.4*(f-f_r)/f_r);
end