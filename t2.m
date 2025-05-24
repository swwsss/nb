% 问题2：霄驰V1机场-千岛湖直飞可行性分析（含最大行李假设）
clc; clear;

% 1. 参数定义（可根据题意适当调整）
g = 9.79;             % 重力加速度 (m/s^2)
rho = 1.18;           % 空气密度 (kg/m^3)
S = 2;                % 迎风面积 (m^2)
Cd = 0.3;             % 阻力系数
eta = 0.8;            % 能效
bat_mass = 200;       % 电池质量 (kg)
bat_den = 300;        % 电池能量密度 (Wh/kg)
E_bat_Wh = bat_mass * bat_den;     % 电池总能量 (Wh)
E_bat_J = E_bat_Wh * 3600;         % 电池总能量 (J)
h = 200;              % 爬升/巡航高度 (m) 可调整
v_climb = 5;          % 垂直爬升/降落速度 (m/s)

% 载重参数（司机+3乘客+4件最大行李，含电池）
m_empty = 800;             % 空车质量 (kg)，含电池不含司机与乘客和行李
m_driver = 65;             % 司机 (kg)
m_passenger = 3*65;        % 3名乘客 (kg)
m_luggage = 4*20;          % 4人各带1件最大行李 (kg)
m = m_empty + m_driver + m_passenger + m_luggage; % 总质量 (kg)

% 2. 经纬度（机场、千岛湖广场，单位：度）
% 萧山国际机场(120.432, 30.229), 千岛湖广场(119.036, 29.604)
lon1 = 120.432; lat1 = 30.229;
lon2 = 119.036; lat2 = 29.604;
R = 6371; % 地球半径 (km)
dlat = deg2rad(lat2 - lat1);
dlon = deg2rad(lon2 - lon1);
a = sin(dlat/2)^2 + cos(deg2rad(lat1))*cos(deg2rad(lat2))*sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
L_actual = R * c * 1000; % 实际直线距离 (米)
L_actual_km = L_actual / 1000;

% 3. 最大航程遍历（与第一题一致，自动找最优速度）
V_kmh = 120:1:160;   % 巡航速度区间 (km/h)
V = V_kmh/3.6;       % 巡航速度 (m/s)
E_up = m * g * h / eta;
E_down = 0.4 * E_up;
E_updown = E_up + E_down;
L_max = zeros(size(V)); % 最大航程 (米)
for i = 1:length(V)
    denom = 0.5 * rho * S * Cd * V(i)^2 / eta;
    L_max(i) = (E_bat_J - E_updown) / denom;
end
L_max_km = L_max / 1000;

% 4. 判断能否直飞（找最大航程覆盖实际距离的最优速度）
idx = find(L_max >= L_actual, 1, "last");
if isempty(idx)
    fprintf('\n--- 霄驰V1无法满载直飞机场-千岛湖广场（安全最大行李方案） ---\n');
    [maxL, maxidx] = max(L_max_km);
    fprintf('最大航程 %.2f km，实际需求 %.2f km\n', maxL, L_actual_km);
    fprintf('最优速度：%.2f km/h\n', V_kmh(maxidx));
    fprintf('建议：需要中途充电或分段飞行。\n');
else
    best_V = V_kmh(idx);
    best_V_ms = V(idx);
    t_cruise = L_actual / best_V_ms;
    t_up = h / v_climb;
    t_down = h / v_climb;
    t_total = (t_cruise + t_up + t_down) / 3600; % 小时

    fprintf('\n--- 霄驰V1可满载直飞机场-千岛湖广场（安全最大行李方案） ---\n');
    fprintf('实际距离：%.2f km\n', L_actual_km);
    fprintf('最优巡航速度：%.2f km/h (%.2f m/s)\n', best_V, best_V_ms);
    fprintf('推荐飞行高度：%.0f m\n', h);
    fprintf('总飞行时间：%.2f 小时\n', t_total);
    fprintf('最大可飞航程：%.2f km\n', L_max(idx)/1000);
end

% 5. 可视化（选做）
% 1. 推荐速度索引
nV = length(V);
energy_used = zeros(1, nV);
energy_remain = zeros(1, nV);
time_total = zeros(1, nV);

for i = 1:nV
    % 能耗模型
    E_cruise = 0.5 * rho * S * Cd * V(i)^2 * L_actual / eta;
    E_total = E_updown + E_cruise;
    energy_used(i) = E_total / 3600 / 1000; % kWh
    energy_remain(i) = E_bat_Wh - E_total / 3600; % Wh
    % 时间模型
    t_cruise = L_actual / V(i);
    t_total = (t_cruise + t_up + t_down) / 3600; % 小时
    time_total(i) = t_total;
end

idx_fast = find(L_max >= L_actual, 1, 'last');   % 最快（最大可行速度）
idx_safe = find(L_max >= L_actual, 1, 'first');  % 最省电（最小可行速度）

% 2. 剩余电量 vs 巡航速度
figure;
plot(V_kmh, energy_remain/1000, 'g-', 'LineWidth',2); grid on; hold on;
if ~isempty(idx_fast)
    plot(V_kmh(idx_fast), energy_remain(idx_fast)/1000, 'ro', 'MarkerSize',10,'LineWidth',2);
end
if ~isempty(idx_safe)
    plot(V_kmh(idx_safe), energy_remain(idx_safe)/1000, 'bo', 'MarkerSize',10,'LineWidth',2);
end
xlabel('巡航速度 (km/h)');
ylabel('剩余电量 (kWh)');
title('飞完全程后剩余电量 vs 巡航速度');
legend('剩余电量','最快安全速度','最省电安全速度');

% 3. 总飞行时间 vs 巡航速度
figure;
plot(V_kmh, time_total, 'm-', 'LineWidth',2); grid on; hold on;
if ~isempty(idx_fast)
    plot(V_kmh(idx_fast), time_total(idx_fast), 'ro', 'MarkerSize',10,'LineWidth',2);
end
if ~isempty(idx_safe)
    plot(V_kmh(idx_safe), time_total(idx_safe), 'bo', 'MarkerSize',10,'LineWidth',2);
end
xlabel('巡航速度 (km/h)');
ylabel('总飞行时间 (小时)');
title('机场-千岛湖总飞行时间 vs 巡航速度');
legend('总飞行时间','最快安全速度','最省电安全速度');

% 4. 剩余电量 vs 总飞行时间（多目标权衡散点图）
figure;
scatter(time_total, energy_remain/1000, 40, V_kmh, 'filled');
colorbar; colormap jet;
xlabel('总飞行时间 (小时)');
ylabel('剩余电量 (kWh)');
title('速度-能耗-时间多目标权衡图');
hold on;
if ~isempty(idx_fast)
    scatter(time_total(idx_fast), energy_remain(idx_fast)/1000, 80, 'r', 'filled');
    text(time_total(idx_fast), energy_remain(idx_fast)/1000, sprintf('最快安全速度%.0fkm/h', V_kmh(idx_fast)), ...
        'VerticalAlignment','bottom','HorizontalAlignment','right');
end
if ~isempty(idx_safe)
    scatter(time_total(idx_safe), energy_remain(idx_safe)/1000, 80, 'b', 'filled');
    text(time_total(idx_safe), energy_remain(idx_safe)/1000, sprintf('最省电速度%.0fkm/h', V_kmh(idx_safe)), ...
        'VerticalAlignment','top','HorizontalAlignment','left');
end
legend('各速度方案','最快安全速度','最省电速度');

% 5. 能耗结构饼图（以最快安全速度为例）
if ~isempty(idx_fast)
    best_V_ms = V(idx_fast);
    E_cruise = 0.5 * rho * S * Cd * best_V_ms^2 * L_actual / eta;
    E_up_pie = m * g * h / eta;
    E_down_pie = 0.4 * E_up_pie;
    pie_data = [E_up_pie, E_cruise, E_down_pie];
    pie_labels = {'起飞爬升','巡航','降落'};
    figure; 
    pie(pie_data);
    legend(pie_labels, 'Location', 'bestoutside');
    title('飞行任务能耗结构（最快安全速度）');
end