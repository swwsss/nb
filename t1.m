% 霄驰V1 eVTOL 最大航程计算脚本
% 适合MATLAB任意常用版本，直接可运行

clc; clear;

%------ 参数设置 ---------
m = 1300;                   % 总质量 (kg)
g = 9.79;                   % 重力加速度 (m/s^2)
rho = 1.18;                 % 空气密度 (kg/m^3)
S = 2;                      % 迎风面积 (m^2)
Cd = 0.3;                   % 阻力系数
eta = 0.8;                  % 能效
bat_mass = 200;             % 电池质量(kg)
bat_den = 300;              % 电池能量密度(Wh/kg)
E_bat_Wh = bat_mass * bat_den;     % 总能量(Wh)
E_bat_J = E_bat_Wh * 3600;         % 总能量(J)
h = 200;                    % 爬升高度(可调，单位m)
v_climb = 5;                % 爬升/降落速度(m/s)
V_kmh = 120:1:160;          % 巡航速度区间(km/h)
V = V_kmh / 3.6;            % 巡航速度(m/s)

%------ 1. 起飞/降落能耗（J） -------
E_up = m * g * h / eta;
E_down = 0.4 * E_up;
E_updown = E_up + E_down;

%------ 2. 计算不同速度下最大航程 -----
L = zeros(size(V));  % 初始化最大航程数组（m）
for i = 1:length(V)
    denom = 0.5 * rho * S * Cd * V(i)^2 / eta;
    L(i) = (E_bat_J - E_updown) / denom;   % 单位：米
end

L_km = L / 1000;     % 单位换算km

%------ 3. 找到最大航程对应速度 ------
[max_L, idx] = max(L_km);
best_V = V_kmh(idx); % 最优速度(km/h)
best_V_ms = V(idx);  % 最优速度(m/s)
best_L = max_L;      % 最优航程(km)

%------ 4. 计算总时间 ------
t_cruise = (best_L * 1000) / best_V_ms;    % 巡航时间(s)
t_up = h / v_climb;                        % 起飞时间(s)
t_down = h / v_climb;                      % 降落时间(s)
t_total = (t_cruise + t_up + t_down) / 3600; % 小时

%------ 5. 输出结果 -------
fprintf('------ 霄驰V1 满载最大航程计算结果 ------\n');
fprintf('总质量：%.0f kg，巡航高度：%.0f m\n', m, h);
fprintf('最优巡航速度：%.2f km/h (%.2f m/s)\n', best_V, best_V_ms);
fprintf('最大航程：%.2f km\n', best_L);
fprintf('总飞行时间：%.2f 小时\n', t_total);
fprintf('起飞+降落总能耗：%.2f kWh (%.2f%%电池)\n', E_updown/3600/1000, E_updown/E_bat_J*100);

%------ 6. 可视化 ------
figure;
plot(V_kmh, L_km, 'b-', 'LineWidth',2); grid on;
xlabel('巡航速度 (km/h)');
ylabel('最大航程 (km)');
title('霄驰V1最大航程 vs 巡航速度');
hold on;
plot(best_V, best_L, 'ro', 'MarkerSize',10,'LineWidth',2);
legend('最大航程','最优速度点');

%------ 7. 能耗结构饼图 ------
E_cruise = 0.5 * rho * S * Cd * best_V_ms^2 * (best_L*1000) / eta; % 巡航能耗(J)
pie_data = [E_up, E_cruise, E_down];
pie_labels = {'起飞爬升','巡航','降落'};

figure;
pie(pie_data);
legend(pie_labels, 'Location', 'bestoutside');
title('霄驰V1飞行任务能耗结构');
set(findobj(gca,'Type','text'),'String','');
percent = pie_data/sum(pie_data)*100;
str = sprintf('起飞爬升: %.1f%%   巡航: %.1f%%   降落: %.1f%%', percent);
annotation('textbox',[0.15 0.05 0.8 0.1],...
    'String',str,'EdgeColor','none','HorizontalAlignment','center','FontSize',12);
