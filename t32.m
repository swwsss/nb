% 问题二：带噪声约束，自动遍历最小可行高度，所有可视化严格仅用噪声达标速度集合
clc; clear;

% 参数定义
g = 9.79;             
rho = 1.18;           
S = 2;                
Cd = 0.3;             
eta = 0.8;            
bat_mass = 200;       
bat_den = 300;        
E_bat_Wh = bat_mass * bat_den;    
E_bat_J = E_bat_Wh * 3600;        
v_climb = 5;          
noise_limit = 55;     % 地面噪声限制

% 载重参数
m_empty = 800;             
m_driver = 65;
m_passenger = 3*65;
m_luggage = 4*20;
m = m_empty + m_driver + m_passenger + m_luggage;

% 经纬度及实际距离
lon1 = 120.432; lat1 = 30.229;
lon2 = 119.036; lat2 = 29.604;
R = 6371; 
dlat = deg2rad(lat2 - lat1);
dlon = deg2rad(lon2 - lon1);
a = sin(dlat/2)^2 + cos(deg2rad(lat1))*cos(deg2rad(lat2))*sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
L_actual = R * c * 1000; % 米
L_actual_km = L_actual / 1000;

% 读取数据并英文化列名
filename = '附件3 巡航速度高度旋翼转速噪声统计数据.xlsx';
data = readtable(filename);
data.Properties.VariableNames = {'xunhangsudu', 'feixinggaodu', 'xuanyizhuansu', 'dimianzaosheng'};

can_direct_fly = false;
heights_all = unique(data.feixinggaodu);

for h = heights_all'
    idx_h = data.feixinggaodu == h & data.dimianzaosheng <= noise_limit;
    V_kmh_plot = data.xunhangsudu(idx_h); % 只筛选噪声合规速度
    dimianzaosheng_vec_plot = data.dimianzaosheng(idx_h);
    if isempty(V_kmh_plot)
        continue;
    end
    V_plot = V_kmh_plot / 3.6;
    E_up = m * g * h / eta;
    E_down = 0.4 * E_up;
    E_updown_plot = E_up + E_down;
    t_up_plot = h / v_climb;
    t_down_plot = h / v_climb;
    nV = length(V_plot);
    L_max_plot = zeros(1, nV);
    for i = 1:nV
        denom = 0.5 * rho * S * Cd * V_plot(i)^2 / eta;
        L_max_plot(i) = (E_bat_J - E_updown_plot) / denom; % m
    end
    idx = find(L_max_plot >= L_actual, 1, "last");
    if ~isempty(idx)
        can_direct_fly = true;
        best_h = h;
        best_V = V_kmh_plot(idx);
        best_V_ms = V_plot(idx);
        best_idx = idx;
        best_noise = dimianzaosheng_vec_plot(idx);
        break;
    end
end

if ~can_direct_fly
    disp('所有高度下都无法满足噪声+航程双重约束，建议采用分段/中途充电方案。');
else
    % 能耗、时间等分析
    t_cruise = L_actual / best_V_ms;
    t_total = (t_cruise + t_up_plot + t_down_plot) / 3600; % 小时
    fprintf('\n--- 霄驰V1可满载直飞机场-千岛湖广场（地面噪声限制） ---\n');
    fprintf('最小可行巡航高度：%.0f m\n', best_h);
    fprintf('最快可行速度（最大速度且航程覆盖）：%.2f km/h (%.2f m/s)\n', best_V, best_V_ms);
    fprintf('对应最大航程：%.2f km\n', L_max_plot(best_idx)/1000);
    fprintf('实际航段距离：%.2f km\n', L_actual_km);
    fprintf('总飞行时间：%.2f 小时\n', t_total);
    fprintf('此速度地面噪声：%.1f dB(A)\n', best_noise);

    % 可视化分析（只用噪声达标速度区间）
    nV = length(V_plot);
    energy_used = zeros(1, nV);
    energy_remain = zeros(1, nV);
    time_total = zeros(1, nV);
    for i = 1:nV
        E_cruise = 0.5 * rho * S * Cd * V_plot(i)^2 * L_actual / eta;
        E_total = E_updown_plot + E_cruise;
        energy_used(i) = E_total / 3600 / 1000; % kWh
        energy_remain(i) = E_bat_Wh - E_total / 3600; % Wh
        t_cruise_i = L_actual / V_plot(i);
        t_total_i = (t_cruise_i + t_up_plot + t_down_plot) / 3600; % 小时
        time_total(i) = t_total_i;
    end

    % “最省电安全速度”：能飞完全程的最小速度
    idx_safe = find(L_max_plot >= L_actual, 1, 'first');
    % “最快安全速度”：能飞完全程的最大速度
    idx_fast = find(L_max_plot >= L_actual, 1, 'last');

    % 剩余电量 vs 巡航速度
    figure;
    plot(V_kmh_plot, energy_remain/1000, 'g-', 'LineWidth',2); grid on; hold on;
    if ~isempty(idx_fast)
        plot(V_kmh_plot(idx_fast), energy_remain(idx_fast)/1000, 'ro', 'MarkerSize',10,'LineWidth',2);
    end
    if ~isempty(idx_safe)
        plot(V_kmh_plot(idx_safe), energy_remain(idx_safe)/1000, 'bo', 'MarkerSize',10,'LineWidth',2);
    end
    xlabel('巡航速度 (km/h)');
    ylabel('剩余电量 (kWh)');
    title(['飞完全程后剩余电量 vs 巡航速度（', num2str(best_h), 'm，噪声约束）']);
    legend('剩余电量','最快安全速度','最省电安全速度');

    % 总飞行时间 vs 巡航速度
    figure;
    plot(V_kmh_plot, time_total, 'm-', 'LineWidth',2); grid on; hold on;
    if ~isempty(idx_fast)
        plot(V_kmh_plot(idx_fast), time_total(idx_fast), 'ro', 'MarkerSize',10,'LineWidth',2);
    end
    if ~isempty(idx_safe)
        plot(V_kmh_plot(idx_safe), time_total(idx_safe), 'bo', 'MarkerSize',10,'LineWidth',2);
    end
    xlabel('巡航速度 (km/h)');
    ylabel('总飞行时间 (小时)');
    title(['机场-千岛湖总飞行时间 vs 巡航速度（', num2str(best_h), 'm，噪声约束）']);
    legend('总飞行时间','最快安全速度','最省电安全速度');

    % 剩余电量 vs 总飞行时间
    figure;
    scatter(time_total, energy_remain/1000, 80, V_kmh_plot, 'filled');
    colorbar; colormap jet;
    xlabel('总飞行时间 (小时)');
    ylabel('剩余电量 (kWh)');
    title(['速度-能耗-时间多目标权衡图（', num2str(best_h), 'm，噪声约束）']);
    hold on;
    if ~isempty(idx_fast)
        scatter(time_total(idx_fast), energy_remain(idx_fast)/1000, 120, 'r', 'filled');
        text(time_total(idx_fast), energy_remain(idx_fast)/1000, ...
            sprintf('最快安全速度%.0fkm/h', V_kmh_plot(idx_fast)), ...
            'VerticalAlignment','bottom','HorizontalAlignment','right');
    end
    if ~isempty(idx_safe)
        scatter(time_total(idx_safe), energy_remain(idx_safe)/1000, 120, 'b', 'filled');
        text(time_total(idx_safe), energy_remain(idx_safe)/1000, ...
            sprintf('最省电速度%.0fkm/h', V_kmh_plot(idx_safe)), ...
            'VerticalAlignment','top','HorizontalAlignment','left');
    end
    legend('各速度方案','最快安全速度','最省电安全速度');

    % 能耗结构饼图（以最快安全速度为例）
    if ~isempty(idx_fast)
        best_V_ms_plot = V_plot(idx_fast);
        E_cruise = 0.5 * rho * S * Cd * best_V_ms_plot^2 * L_actual / eta;
        E_up_pie = m * g * best_h / eta;
        E_down_pie = 0.4 * E_up_pie;
        pie_data = [E_up_pie, E_cruise, E_down_pie];
        pie_labels = {'起飞爬升','巡航','降落'};
        figure; 
        pie(pie_data);
        legend(pie_labels, 'Location', 'bestoutside');
        title(['飞行任务能耗结构（最快安全速度，', num2str(best_h), 'm，噪声约束）']);
    end
end