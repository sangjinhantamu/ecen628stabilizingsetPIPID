function [lower_bound, upper_bound]=sangjin_stabilizing_set_rational_pi_01(Kp, p_1_w_poly, p_2_w_poly, q_1_w_poly, q_2_w_poly, sig_n, sig_m, sig_nu, Ki_min, Ki_max)
%Kp = 1;
syms w 'positive';
p_1_w = poly2sym(p_1_w_poly, w);
p_2_w = poly2sym(p_2_w_poly, w);
q_1_w = poly2sym(q_1_w_poly, w);
q_2_w = poly2sym(q_2_w_poly, w);
q_w = q_1_w + Kp*q_2_w;
w_q_Kp = roots(sym2poly(q_w));
w_q_Kp(any((w_q_Kp < 0),2))=[];
w_q_Kp(any(imag(w_q_Kp),2))=[];
% test
% w_q_Kp = [w_q_Kp; w_q_Kp(end)]
dup_w_q_Kp = sum(w_q_Kp==w_q_Kp');
w_q_Kp(any((mod(dup_w_q_Kp,2)'==0),2))=[];
w_q_Kp = sort(w_q_Kp);
if(mod(sig_n+sig_m,2)~=0)
    w_q_Kp = [w_q_Kp; inf];
end
% the book says sgn[q(0,-18)] = -1 but it should be q(0+,-18)) and the
% smallest positive number in MATLAB is 'realmin'.
sgn_j = sign(subs(q_w,w,realmin));
str_comb = 2*ones(1,numel(w_q_Kp));
str_base = [-1 1];
comb_i = str_base(fullfact(str_comb));
str_i = ones(1,numel(w_q_Kp));
for idx = 1:numel(w_q_Kp)
    if(isfinite(w_q_Kp(idx)) && w_q_Kp(idx) >0)
        isfinitepos = 2;
    else
        isfinitepos = 1;
    end
    str_i(idx) = (-1)^(idx+1) * isfinitepos;
end
comb_i_temp=comb_i.*str_i;
admissible_strs = comb_i(sgn_j*sum(comb_i_temp,2) == sig_nu,:);
num_strs = size(admissible_strs,1);
num_ineqs = size(admissible_strs,2);

if(num_strs >1)
    if Kp>0
        sangjin=0;
    end
    %error('Ki intervals are more than one.');
    lower_bound=NaN(num_strs,1);
    upper_bound=NaN(num_strs,1);
    for idx=1:num_strs
        % Ax > b
        A_linineq=zeros(num_ineqs,1);
        b_linineq=zeros(num_ineqs,1);
        for idy=1:num_ineqs
            w_freq = w_q_Kp(idy);
            p_1_w_freq = double(subs(p_1_w,w,w_freq));
            p_2_w_freq = double(subs(p_2_w,w,w_freq));
            if(admissible_strs(idx,idy)*p_2_w_freq < 0)
                A_linineq(idy,1)=-1;
                b_linineq(idy)=((admissible_strs(idx,idy))*p_1_w_freq)/((admissible_strs(idx,idy))*p_2_w_freq);
            else
                A_linineq(idy,1)=1;
                b_linineq(idy)=((admissible_strs(idx,idy))*p_1_w_freq)/((((-1)*admissible_strs(idx,idy)))*p_2_w_freq);
            end
        end
        Ab_linineq=[A_linineq, b_linineq];
        A_linineq(any(isnan(Ab_linineq),2),:)=[];
        b_linineq(any(isnan(Ab_linineq),2),:)=[];
        num_bounds = length(A_linineq);
        upper_bounds = NaN(num_bounds,1);
        lower_bounds = NaN(num_bounds,1);
        for idx_tmp=1:num_bounds
            if(A_linineq(idx_tmp,1) > 0)
                lower_bounds(idx_tmp,1) = b_linineq(idx_tmp,1);
            else
                upper_bounds(idx_tmp,1) = -b_linineq(idx_tmp,1);
            end
        end
        lower_bounds(any(isnan(lower_bounds),2),:)=[];
        upper_bounds(any(isnan(upper_bounds),2),:)=[];
        if(isempty(lower_bounds))
            lower_bound(idx,1)=Ki_min;
        else
            lower_bound(idx,1) = max(lower_bounds); %min(lower_bounds);
        end
        if(isempty(upper_bounds))
            upper_bound(idx,1)=Ki_max;
        else
            upper_bound(idx,1) = min(upper_bounds); %max(upper_bounds);
        end
        if(lower_bound(idx,1) >= upper_bound(idx,1))
            lower_bound(idx,1)=NaN;
            upper_bound(idx,1)=NaN;
        end
    end
    lower_bound(any(isnan(lower_bound),2),:)=[];
    upper_bound(any(isnan(upper_bound),2),:)=[];
    if(isempty(lower_bound))
        lower_bound=NaN;
    end
    if(isempty(upper_bound))
        upper_bound=NaN;
    end
elseif(isempty(admissible_strs))
    warning('No Ki interval');
    upper_bound = NaN;
    lower_bound = NaN;
    return
else
    if Kp>0
        sangjin=0;
    end
    idx=1;
    % Ax > b
    A_linineq=zeros(num_ineqs,1);
    b_linineq=zeros(num_ineqs,1);
    for idy=1:num_ineqs
        w_freq = w_q_Kp(idy);
        p_1_w_freq = double(subs(p_1_w,w,w_freq));
        p_2_w_freq = double(subs(p_2_w,w,w_freq));
        if(admissible_strs(idx,idy)*p_2_w_freq < 0)
            A_linineq(idy,1)=-1;
            b_linineq(idy)=((admissible_strs(idx,idy))*p_1_w_freq)/((admissible_strs(idx,idy))*p_2_w_freq);
        else
            A_linineq(idy,1)=1;
            b_linineq(idy)=((admissible_strs(idx,idy))*p_1_w_freq)/((((-1)*admissible_strs(idx,idy)))*p_2_w_freq);
        end
    end
    Ab_linineq=[A_linineq, b_linineq];
    A_linineq(any(isnan(Ab_linineq),2),:)=[];
    b_linineq(any(isnan(Ab_linineq),2),:)=[];
    num_bounds = length(A_linineq);
    upper_bounds = NaN(num_bounds,1);
    lower_bounds = NaN(num_bounds,1);
    for idx=1:num_bounds
        if(A_linineq(idx,1) > 0)
            lower_bounds(idx,1) = b_linineq(idx,1);
        else
            upper_bounds(idx,1) = -b_linineq(idx,1);
        end
    end
    lower_bounds(any(isnan(lower_bounds),2),:)=[];
    upper_bounds(any(isnan(upper_bounds),2),:)=[];
    if(isempty(lower_bounds))
        lower_bound=Ki_min;
    else
        lower_bound = max(lower_bounds); %min(lower_bounds);
    end
    if(isempty(upper_bounds))
        upper_bound=Ki_max;
    else
        upper_bound = min(upper_bounds); %max(upper_bounds);
    end
    if(lower_bound >= upper_bound)
        lower_bound=NaN;
        upper_bound=NaN;
    end
end
