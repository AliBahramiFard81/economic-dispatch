%% by Ali Bahrami Fard
% alibt1313@gmail.com
% ECONOMIC DISPATCH
clc;
clear;
close all;
format short g;

%% Inputs
syms p1 p2 p3 L; 

H = [510 + 7.2 * p1 + 0.00142 * p1 ^ 2;
     310 + 7.85 * p2 + 0.00194 * p2 ^ 2;
     78 + 7.97 * p3 + 0.00482 * p3 ^ 2];

H_cost = [1.1;
           1;
           1;];

power_limits = [150,600;
                100,400;
                50,200;];

pload = 850;

vars = symvar(H);

%% Calculating Functions
F = H;
for i=1:length(H)
   F(i,1) = H(i,1) * H_cost(i,1);
end

%% Calculating Lambdas
F_diff = F;
for i=1:length(F_diff)
   F_diff(i,1) = diff(F(i,1)) == L;
end
F_diff(length(F_diff) + 1,1) = sum(vars) == pload;

result = 0;
result = 0;
results = 0;


%% Checking if the values are in range
% the F that is out of limits       new power value     0 == min, 1 == max,
% in range == 2
check_limits = [];
for i=2:length(results)
    if results(i,1) < power_limits(i - 1,1)
        check_limits(i - 1,1) = i -1;
        check_limits(i - 1,2) = power_limits(i - 1,1);
        check_limits(i - 1,3) = 0;
    elseif results(i,1) > power_limits(i - 1,2)
        check_limits(i - 1,1) = i -1;
        check_limits(i - 1,2) = power_limits(i - 1,2);
        check_limits(i - 1,3) = 1;
    end
end

for i=1:length(check_limits)
    if check_limits(i,1) == 0
        check_limits(i,1) = i;
        check_limits(i,2) = pload - sum(check_limits(:,2));
        check_limits(i,3) = 2;
    end
end

%% Correcting the out of limit plants (calculating new lambdas)
new_lambdas = zeros(length(check_limits),1);
if ~isempty(check_limits)
    for i=1:length(check_limits)
        new_lambdas(i,1) = double(solve(subs(F_diff(check_limits(i,1),1),check_limits(i,2)),L));
    end
end

%% Checking if the lambdas are in range
out_range_lambda_index = [];
in_range_lambda_index = [];
if ~isempty(new_lambdas)
    in_range_lambda_index = find(check_limits(:,3) > 1);
    in_range_lambda = new_lambdas(in_range_lambda_index,1);
    index = 0;
    for i=1:length(new_lambdas)
        if check_limits(i,3) == 0
            if new_lambdas(i,1) < in_range_lambda
                index = index + 1;
                out_range_lambda_index(index,1) = i;
            end
        end
        if check_limits(i,3) == 1
            if new_lambdas(i,1) > in_range_lambda
                index = index + 1;
                out_range_lambda_index(index,1) = i;
            end
        end
    end
end

%% Solving ED for those that aren't in range
new_fdiff = F_diff;
remaining_power = [];
if ~isempty(out_range_lambda_index)
    in_out_range_syms = vars;
    for i=1:length(out_range_lambda_index)
        new_fdiff(length(new_fdiff) + i,1) = F_diff(out_range_lambda_index(i),1);
        in_out_range_syms(1,length(in_out_range_syms) + 1) = vars(1,out_range_lambda_index(i));
    end
    for i=1:length(in_range_lambda_index)
        new_fdiff(length(new_fdiff) + i,1) = F_diff(in_range_lambda_index(i),1);
        in_out_range_syms(1,length(in_out_range_syms) + 1) = vars(1,in_range_lambda_index(i));
    end
    in_out_range_syms = in_out_range_syms(1,length(vars) + 1:length(in_out_range_syms));
    new_fdiff(length(new_fdiff) + 1,1) = sum(in_out_range_syms) == sum(check_limits(out_range_lambda_index,2)) + sum(check_limits(in_range_lambda_index,2));
    new_fdiff = new_fdiff(length(F_diff) + 1:length(new_fdiff),1);

    result = vpasolve(new_fdiff);
    result = struct2array(result);
    results = zeros(length(vars),1);
    for i=1:length(result)
        results(i,1) = double(result(i));
    end
    disp(results(1,1));
    disp(transpose(results(2:end,1)))
end

if isempty(check_limits)
    result = vpasolve(F_diff);
    result = struct2array(result);
    results = zeros(length(vars),1);
    for i=1:length(result)
        results(i,1) = double(result(i));
    end
    disp(results(1,1));
    disp(transpose(results(2:end,1)))
end

%% Calculating Ft
Ft = 0;
for i=2:length(results)
    Ft = Ft + subs(F(i-1,1), results(i));
end
disp(double(Ft))
disp('----------------------------------------------------')

%% Piecewise linear const function
number_of_pieces = 11;
pieced_power_limit = [];
index = 1;
for i=1:length(power_limits)
    power_limit = power_limits(i,:);
    power_difference = (power_limit(1,2) - power_limit(1,1)) / number_of_pieces;
    prev_differencce = power_limit(1,1);
    for x=1:number_of_pieces + 1
        pieced_power_limit(index,x) = prev_differencce;
        prev_differencce = prev_differencce + power_difference;
    end
    index = index + 1;
end

%% Calculating F
F_piece = [];
for i=1:length(pieced_power_limit(:,1))
    for x=1:length(pieced_power_limit(1,:))
       F_piece(i,x) = subs(F(i,1),pieced_power_limit(i,x));
    end
end

%% Calculating S
S = [];
for i=1:length(pieced_power_limit(:,1))
    power_limit = power_limits(i,:);
    power_difference = (power_limit(1,2) - power_limit(1,1)) / number_of_pieces;
    for x=1:length(pieced_power_limit(1,:)) - 1
        S(i,x) = (F_piece(i,x + 1) - F_piece(i,x)) / power_difference;
    end
end
disp(pieced_power_limit)
disp(S)

%% Calculating F for each piece
F_x = [];
for i=1:length(pieced_power_limit(:,1))
    for x=1:length(pieced_power_limit(i,:))
        F_x(i,x) = double(subs(F(i,1),pieced_power_limit(i,x)));
    end
end
disp(F_x)

%% Ploting the pieces
figure(1);
plot(pieced_power_limit(1,:), F_x(1,:),'DisplayName','Y -X','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','*')
xlabel("P")
ylabel("F")
yticks(F_x(1,:))
xticks(pieced_power_limit(1,:))
title("Unit 1")

figure(2);
plot(pieced_power_limit(2,:), F_x(2,:),'DisplayName','Y -X','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','*')
xlabel("P")
ylabel("F")
yticks(F_x(2,:))
xticks(pieced_power_limit(2,:))
title("Unit 2")

figure(3);
plot(pieced_power_limit(3,:), F_x(3,:),'DisplayName','Y -X','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','*')
xlabel("P")
ylabel("F")
yticks(F_x(3,:))
xticks(pieced_power_limit(3,:))
title("Unit 3")
%% Calculating power for each plant
s_rows = [];
index = 1;
sorted_s = sort(S(:));
smallest_s = min(sorted_s);
calculated_p = zeros(1,length(H));

for i=1:length(sorted_s)
    [r,c] = find(S == sorted_s(i));
    s_rows(i,1) = r; 
end
for i=1:length(s_rows)
    power_limit = power_limits(s_rows(i),:);
    power_difference = (power_limit(1,2) - power_limit(1,1)) / number_of_pieces;
    prev_differencce = power_limit(1,1);
    if sum(calculated_p) < pload
        if calculated_p(s_rows(i)) == 0
            calculated_p(1,s_rows(i)) =  prev_differencce;
        else
            calculated_p(1,s_rows(i)) = calculated_p(1,s_rows(i)) + power_difference;
        end
    end
    if sum(calculated_p) > pload
        calculated_p(1,s_rows(i)) = calculated_p(1,s_rows(i)) - power_difference;
        calculated_p(1,s_rows(i)) = calculated_p(1,s_rows(i)) + pload - sum(calculated_p);
        break; 
    end
end

if sum(calculated_p) < pload
    [r,c] = find(S == smallest_s);
    calculated_p(1,r) = calculated_p(1,r) + pload - sum(calculated_p);
end
disp(calculated_p)
Ft_piece = 0;
for i=1:length(calculated_p)
    Ft_piece = Ft_piece + subs(F(i,1), calculated_p(i));
end
disp(double(Ft_piece))
disp('----------------------------------------------------')

%% Lambda iteration method
initial_lambda = 8;

% first iterration
F_diff_it = F;
for i=1:length(F_diff_it)
   F_diff_it(i,1) = diff(F(i,1)) == initial_lambda;
end
result_it = vpasolve(F_diff_it);
result_it = struct2array(result_it);
results_it = zeros(length(vars),1);
for i=1:length(result_it)
    results_it(i,1) = double(result_it(i));
end

first_epsilon = sum(results_it) - pload;

% second iterration
second_lambda = 0;
if first_epsilon < 0
    second_lambda = initial_lambda * 1.1;
else 
    second_lambda = initial_lambda * 0.9;
end

disp(first_epsilon)
disp(second_lambda);
for i=1:length(F_diff_it)
   F_diff_it(i,1) = diff(F(i,1)) == second_lambda;
end
result_it = vpasolve(F_diff_it);
result_it = struct2array(result_it);
results_it = zeros(length(vars),1);
for i=1:length(result_it)
    results_it(i,1) = double(result_it(i));
end
second_epsilon = sum(results_it) - pload;

disp(second_epsilon)
counter = 0;
while true
    M = (second_epsilon - first_epsilon) / (second_lambda - initial_lambda);
    B = -1 * (M * second_lambda) + second_epsilon;
    
    initial_lambda = second_lambda;
    second_lambda = (-1 * B) / M;

    F_diff_it = F;
    for i=1:length(F_diff_it)
       F_diff_it(i,1) = diff(F(i,1)) == second_lambda;
    end
    result_it = vpasolve(F_diff_it);
    result_it = struct2array(result_it);
    results_it = zeros(length(vars),1);
    for i=1:length(result_it)
        results_it(i,1) = double(result_it(i));
    end
    
    first_epsilon = second_epsilon;
    second_epsilon = sum(results_it) - pload;
    counter = counter + 1;
    disp(results_it)
    disp(second_epsilon)
    disp(second_lambda)
    disp(counter + 2)
    if second_epsilon < 1
        break;
    end
end

disp('----------------------------------------------------')
%% Composite generation cost

% Calculating lambda min and max
lambdas = [];
F_diff_com = F;
for i=1:length(F_diff_com)
   F_diff_com(i,1) = diff(F(i,1));
end

for i=1:length(F_diff_com)
    for x=1:length(power_limits(1,:))
        lambdas(i,x) = double(subs(F_diff_com(i,1),power_limits(i,x)));
    end
end
number_of_pieces = 10;
lambda_min = min(lambdas(:,1));
lambda_max = max(lambdas(:,2));
delta_lambda = (lambda_max - lambda_min) / number_of_pieces;
lambda_pieces = zeros(1,number_of_pieces);
lambda_pieces(1,1) = lambda_min;

for i=2:number_of_pieces + 1
     lambda_pieces(1,i) = lambda_pieces(1,i - 1) + delta_lambda;
end


ps = [];
fs = [];
for i=1:length(lambda_pieces)
    pi = 0;
    fi = 0;
    for x=1:length(F_diff_com)
        px = double(vpasolve(F_diff_com(x,1) == lambda_pieces(1,i)));
        if px < power_limits(x,1)
            px = power_limits(x,1);
        elseif px > power_limits(x,2)
            px = power_limits(x,2);
        end
        fx = double(subs(F(x,1),px));
        pi = pi + px;
        fi = fi + fx;
    end
    ps(1,i) = pi;
    fs(1,i) = fi;
end
disp(lambda_pieces)
disp(ps)
disp(fs)

disp('----------------------------------------------------')

%%  BASE POINT AND PARTICIPATION FACTORS
part_f = [];
F_diff_bp = F;
new_pload = 900;
delta_pload = new_pload - pload;
for i=1:length(F_diff_bp)
   F_diff_bp(i,1) = diff(diff(F(i,1)));
end

delta_p_d = 0;
for i=1:length(F_diff_bp)
    delta_p_d = delta_p_d + (1 / F_diff_bp(i,1));
end
delta_p_d = double(delta_p_d);

for i=1:length(F_diff_bp)
    delta_p_i = 1 / double(F_diff_bp(i,1));
    part_f(1,i) = delta_p_i / delta_p_d;
end

new_p = [];
for i=2:length(results)
    new_p(1,i - 1) = results(i) + part_f(1,i - 1) * delta_pload;
end
disp(delta_pload)
disp(part_f)
disp(new_p)