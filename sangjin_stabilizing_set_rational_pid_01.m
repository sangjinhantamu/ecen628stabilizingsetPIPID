% Fix Kp = - 18 from the allowable range, get the real, nonnegative,
% distict finite zeros of q(w,-Kp) with odd multiplicities, define the last
% omega to be infty, find admissible string, find the inequalities for the
% stabilizing (k_i, K_d) values.
close all;
% Kp=10;

figure; hold on;
for Kp=-6:0.5:1.5
q_w = q_1_w + Kp*q_2_w;
w_q_Kp = roots(sym2poly(q_w));
w_q_Kp(any(imag(w_q_Kp),2),:)=[]; % real
w_q_Kp(any((w_q_Kp < 0),2),:)=[]; % positive
uniqueVals=unique(w_q_Kp);
valCount=hist(w_q_Kp,uniqueVals)';
uniqueVals(any((mod(valCount,2)==0),2),:)=[]; % odd multiplicity
w_q_Kp = sort(uniqueVals);
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

for idx=1:num_strs
    % Ax > b
    A_linineq=zeros(num_ineqs,2);
    b_linineq=zeros(num_ineqs,1);
    for idy=1:num_ineqs
        w_freq = w_q_Kp(idy);
        p_1_w_freq = double(subs(p_1_w,w,w_freq));
        p_2_w_freq = double(subs(p_2_w,w,w_freq));
        if(admissible_strs(idx,idy)*p_2_w_freq < 0)
            A_linineq(idy,1)=-1;
            A_linineq(idy,2)=w_freq^2;
            b_linineq(idy)=((admissible_strs(idx,idy))*p_1_w_freq)/((admissible_strs(idx,idy))*p_2_w_freq);
        else
            A_linineq(idy,1)=1;
            A_linineq(idy,2)=-w_freq^2;
            b_linineq(idy)=((admissible_strs(idx,idy))*p_1_w_freq)/((((-1)*admissible_strs(idx,idy)))*p_2_w_freq);
        end
    end
    admissible_strs(idx,:)
    Ab_linineq=[A_linineq, b_linineq]
    A_linineq(any(isnan(Ab_linineq),2),:)=[];
    b_linineq(any(isnan(Ab_linineq),2),:)=[];
    A_linineq=[A_linineq, ones(size(A_linineq,1),1)];
    b_linineq=b_linineq+Kp*ones(size(b_linineq,1),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% careful tuning is required
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% here!!!
    plotregion_sang(A_linineq,b_linineq,[-200 -200 Kp],[200 200 Kp],[1.0 1.0 1.0]);
%     h = gcf; %current figure handle
%     axesObjs = get(h, 'Children');
%     dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
%     objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
%     xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
%     ydata = get(dataObjs, 'YData');
    % zdata = get(dataObjs, 'ZData');
end

end
axis([-11 0 -11 5 -6 1.5]);
% axis([-10 0 12.95 13.05 9 18]);
%  xlabel('K_i');
%  ylabel('K_d');
%  zlabel('K_p');
 grid on;
% x_min=min(cell2mat(xdata));
% x_max=max(cell2mat(xdata));
% y_min=min(cell2mat(ydata));
% y_max=max(cell2mat(ydata));
% axis([x_min x_max y_min y_max]);
