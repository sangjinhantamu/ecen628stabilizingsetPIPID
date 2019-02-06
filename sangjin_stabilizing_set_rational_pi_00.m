% Compute stabilizing set of PI controllers for finite dimensional LTI
% rational transfer function using signature method.
clear; close all;
% filename = 'stab_pi_plant_00.mat';
% N = [1 -2];         % Numerador
% D = [1 4 3];        % Denominador
% % 'stab_pi_plant_01.mat'
% N = 2;
% D = [1 0 2];
% % 'stab_pi_plant_02.mat'
% N = [1 0.2 2.26];
% D = [1 0.9 2.49 1.145];
% filename = 'stab_pi_plant_03.mat';
% N = [1 -3];         % Numerador
% D = [1 9 -10];        % Denominador
% filename = 'stab_pi_plant_04.mat';
% N = [1 -1];         % Numerador
% D = [1 0.8 -0.2];        % Denominador
% filename = 'stab_pi_plant_05.mat';
% N = [10 9 362.4 36.16];         % Numerador
% D = [2 2.7255 138.4292 156.471 637.6472 360.1779];        % Denominador
% filename = 'stab_pi_plant_06.mat';
% N = [1 0.26 -0.2];         % Numerador s2 + 0.26s - 0.2
% D = [1 0.9 2.49 1.145];        % Denominador s3 + 0.9s2 + 2.49s + 1.145
% filename = 'stab_pi_plant_07.mat';
% N = [1 0.26 0.2];         % Numerador s2 + 0.26s + 0.2
% D = [1.0000    1.7000    5.6200    5.3060    6.9169    2.7595];        % Denominador s3 + 0.9s2 + 2.49s + 1.145
% filename = 'stab_pi_plant_08.mat';
% N = [1 -2];         % Numerador s2 + 0.26s + 0.2
% %D = [1 -16.933 41.799];        % Denominador s3 + 0.9s2 + 2.49s + 1.145
% D = [1 -22.89 59.67];
% filename = 'stab_pi_plant_09.mat';
% N = [1 -1];         % Numerador
% D = [1 2 3];        % Denominador
% filename = 'stab_pi_plant_10.mat';
% N = [1 -1];         % Numerador
% D = [1 2 3];        % Denominador
% filename = 'stab_pi_plant_11.mat';
% N = [1 -0.5];         % Numerador (s-0.5)/(s^2 + s + 2.9)
% D = [1 1 3];        % Denominador
filename = 'stab_pi_plant_12.mat';
N = [-1 1];         % Numerador (s-0.5)/(s^2 + s + 2.9)
D = [1 1 2];        % Denominador


w = sym('w','positive');
s = sym('s');
N_s = poly2sym(N,s);
D_s = poly2sym(D,s);

Ki = sym('Ki', 'real');
Kp = sym('Kp', 'real');

delta_s = s*D_s + Ki*N_s + Kp*s*N_s;
N_ns = subs(N_s, s, -s);

N_flip = fliplr(N);
N_e_flip = N_flip(1,1:2:end);
if(length(N_flip) > 1)
    N_o_flip = N_flip(1,2:2:end);
else
    N_o_flip = 0;
end
N_e = fliplr(N_e_flip);
N_o = fliplr(N_o_flip);
N_e_ss = subs(poly2sym(N_e,s),s,s^2);
N_o_ss = subs(poly2sym(N_o,s),s,s^2);

D_flip = fliplr(D);
D_e_flip = D_flip(1,1:2:end);
if(length(D_flip) > 1)
    D_o_flip = D_flip(1,2:2:end);
else
    D_o_flip = 0;
end
D_e = fliplr(D_e_flip);
D_o = fliplr(D_o_flip);
D_e_ss = subs(poly2sym(D_e,s),s,s^2);
D_o_ss = subs(poly2sym(D_o,s),s,s^2);

nu_jw = subs(delta_s*N_ns,s,1i*w);
nu_jw_real = simplify(real(nu_jw));
nu_jw_real_ki_coeffs = coeffs(nu_jw_real,Ki);
p_1_w = expand(nu_jw_real_ki_coeffs(1));
p_2_w = expand(nu_jw_real_ki_coeffs(end));

nu_jw_imag = simplify(imag(nu_jw));
nu_jw_imag_kp_coeffs = coeffs(nu_jw_imag,Kp);
q_1_w = expand(nu_jw_imag_kp_coeffs(1));
q_2_w = expand(nu_jw_imag_kp_coeffs(2));

root_N_ns = roots(sym2poly(N_ns));
sig_n = length(coeffs(D_s,s,'all')) - 1;
sig_m = length(coeffs(N_s,s,'all')) - 1;
sig_z_pos = sum(real(root_N_ns(:)) < 0);
sig_nu = sig_n - sig_m + 1 + 2*sig_z_pos;

figure;
ezplot(-q_1_w/q_2_w, [0 10]);
h = title('$$K_p = -\frac{q_1(\omega)}{q_2(\omega)}, ~~ \sigma(\nu) = $$','interpreter','latex');
origtitle = get(h,'String');
signature_number = num2str(sig_nu);
set(h,'String',[origtitle ' ' signature_number])
xlabel('\omega (frequency in rad)');
ylabel('K_p');

% 'stab_pi_plant_00.mat'
% Kp_range = [1.5-realmin, 1.49:-0.01:-3.99, -3.995];
% 'stab_pi_plant_01.mat'
% Kp_range = [5-realmin, 4.9:-0.01:-1.09];
% 'stab_pi_plant_02.mat'
% Kp_range = [-0.1:-0.01:-0.5];
% 'stab_pi_plant_03.mat'
% Kp_range = [-3.34-realmin, -3.35:-0.01:-11.41];
% Kp_range = linspace(6.3, -6, 1000);
% 'stab_pi_plant_04.mat'
% Kp_range = linspace(-0.2, -1, 1000);% [-0.2-realmin, -0.21:-0.01:-1.79];
% 'stab_pi_plant_05.mat'
% Kp_range = linspace(-5, 35,1000);
% 'stab_pi_plant_06.mat'
% Kp_range = linspace(-0.64, 5.725,1000);
% 'stab_pi_plant_07.mat'
% Kp_range = linspace(-5, 5,1000);
% 'stab_pi_plant_08.mat'
% Kp_range = linspace(21, 29.8,1000);
% 'stab_pi_plant_09.mat' ??
% Kp_range = linspace(-3, 3,200);
% 'stab_pi_plant_10.mat' ??
% Kp_range = linspace(-3, 3,200);
% 'stab_pi_plant_11.mat' ??
% Kp_range = linspace(-1.5, 6,200);
% 'stab_pi_plant_12.mat'
Kp_range = linspace(-2, 2,200);



Ki_bounds = NaN(numel(Kp_range),2);
if N(end) > 0
    Ki_min = 0;
    Ki_max = 1000;
else
    Ki_min = -1000;
    Ki_max = 0;
end
p_1_w_poly = sym2poly(p_1_w);
p_2_w_poly = sym2poly(p_2_w);
q_1_w_poly = sym2poly(q_1_w);
q_2_w_poly = sym2poly(q_2_w);
Ki_bound = NaN(1,2);
Kp_range_extra = [];
Ki_bounds_extra =[];
for idx=1:numel(Kp_range)
    [lb, ub]=sangjin_stabilizing_set_rational_pi_01(Kp_range(idx), p_1_w_poly, p_2_w_poly, q_1_w_poly, q_2_w_poly, sig_n, sig_m, sig_nu, Ki_min, Ki_max);
    Ki_bound = [lb, ub];
    if(size(lb,1)>1)
        warning('more than two Ki intervals.');
        Ki_bounds(idx,:) = Ki_bound(1,:);
        Ki_bound(1,:)=[];
        for idy=1:size(Ki_bound,1)
            Kp_range_extra = [Kp_range_extra; Kp_range(idx)];
            Ki_bounds_extra = [Ki_bounds_extra; Ki_bound(idy,:)];
        end
    else
        Ki_bounds(idx,:) = Ki_bound;
    end
end
Kp_range=Kp_range';
Kp_range(any(isnan(Ki_bounds),2),:)=[];
Ki_bounds(any(isnan(Ki_bounds),2),:)=[];
figure;
h = fill([Ki_bounds(:,1)' flipud(Ki_bounds(:,2))'],[Kp_range' flipud(Kp_range)'],'y');
xlabel('K_i');
ylabel('K_p');
axis([-2, 0, -2, 6]);
save(filename,'Ki_bounds','Kp_range');
% clear;
% close all;
