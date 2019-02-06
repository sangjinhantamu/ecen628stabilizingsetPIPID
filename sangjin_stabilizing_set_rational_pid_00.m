% Compute stabilizing set for finite dimensional LTI rational transfer
% function.
clear all;
%close all;
% plant 1
%N = [1 -2 -1 -1];
%D = [1 2 32 26 65 -8 1];

% filename = 'stab_pid_plant_02.mat';
% N = [10 9 362.4 36.16]; % Numerador
% D = [2 2.7255 138.4292 156.471 637.6472 360.1779]; % Denominador
%N = [1 0.26 0.2];         % Numerador s2 + 0.26s + 0.2
%D = [1 0.9 2.49 1.145];        % Denominador s3 + 0.9s2 + 2.49s + 1.145
% filename = 'stab_pid_plant_03.mat';
% N = [1 0.26 0.2];         % Numerador s2 + 0.26s + 0.2
% D = [1.0000    1.7000    5.6200    5.3060    6.9169    2.7595];        % Denominador s3 + 0.9s2 + 2.49s + 1.145
filename = 'stab_pid_plant_04.mat';
N = [1 -2];         % Numerador
D = [1 4 3];        % Denominador
% filename = 'stab_pid_plant_05.mat';
% N = [1 3];         % Numerador
% D = [1 9 -10];        % Denominador
% filename = 'stab_pid_plant_06.mat';
% N = [1 -1];         % Numerador
% D = [1 -10 7 18];        % Denominador
% filename = 'stab_pid_plant_07.mat';
% N = [1 -1];         % Numerador
% D = [1 -8 5 14];        % Denominador
% filename = 'stab_pid_plant_08_prelim_01.mat';
% N = [1 -2];         % Numerador
% D = [1 -46.5 86 133.5];        % Denominador


w = sym('w','positive');
s = sym('s');
N_s = poly2sym(N,s);
D_s = poly2sym(D,s);
Ki = sym('Ki', 'real');
Kp = sym('Kp', 'real');
Kd = sym('Kd', 'real');

delta_s = s*D_s + (Ki + Kd*s^2)*N_s + Kp*s*N_s;
N_ns = subs(N_s, s, -s);
nu_jw = subs(delta_s*N_ns,s,1i*w);
nu_jw_real = simplify(real(nu_jw));
nu_jw_real_ki_coeffs = coeffs(nu_jw_real,Ki);
p_2_w = expand(nu_jw_real_ki_coeffs(end));
nu_jw_real_kd_ki_coeffs = coeffs(nu_jw_real_ki_coeffs(1),Kd);
p_1_w = expand(nu_jw_real_kd_ki_coeffs(1));

nu_jw_imag = simplify(imag(nu_jw));
nu_jw_imag_kp_coeffs = coeffs(nu_jw_imag,Kp);
q_1_w = expand(nu_jw_imag_kp_coeffs(1));
q_2_w = expand(nu_jw_imag_kp_coeffs(2));

root_N_ns = roots(sym2poly(N_ns));
sig_n = length(coeffs(D_s,s,'all')) - 1;
sig_m = length(coeffs(N_s,s,'all')) - 1;
sig_z_pos = sum(real(root_N_ns(:)) < 0);
sig_nu = sig_n - sig_m + 1 + 2*sig_z_pos;
ezplot(-q_1_w/q_2_w,[0 10]);

h = title('$$K_p = -\frac{q_1(\omega)}{q_2(\omega)}, ~~ \sigma(\nu) = $$','interpreter','latex');
origtitle = get(h,'String');
signature_number = num2str(sig_nu);
set(h,'String',[origtitle ' ' signature_number])
xlabel('\omega (frequency in rad)');
ylabel('K_p');


