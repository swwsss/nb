% 霄驰V1 eVTOL 最大航程计算脚本（考虑地面噪声限制，变量名英文化，直接可运行）
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
v_climb = 5;                % 爬升/降落速度(m/s)
noise_limit = 55;           % 昼间地面噪声限制

%------ 1. 读取噪声数据 ---------
filename = '附件3 巡航速度高度旋翼转速噪声统计数据.xlsx';
data = readtable(filename);
% 替换变量名为英文
data.Properties.VariableNames = {'xunhangsudu', 'feixinggaodu', 'xuanyizhuansu', 'dimianzaosheng'};

%------ 2. 遍历高度，筛选噪声达标的最大航程 ------
max_L = 0;
best_V = NaN;
best_h = NaN;
best_idx = NaN;
L_record = [];
V_record = [];
H_record = [];
for h = 200:10:500 % 高度步进可改
    % 只考虑该高度、且地面噪声达标的数据    
idx = data.feixinggaodu == h & data.dimianzaosheng <= noise_limit;    
if ~any(idx)        
continue; % 该高度下没有达标速度    
end    
speeds = data.xunhangsudu(idx);    
for v_kmh = speeds'        
V = v_kmh/3.6; % m/s        
        % 起飞/降落能耗
        E_up = m * g * h / eta;
        E_down = 0.4 * E_up;
        E_updown = E_up + E_down;        
        % 最大航程        
denom = 0.5 * rho * S * Cd * V^2 / eta;        
L = (E_bat_J - E_updown) / denom; % m
        L_km = L / 1000;        
        % 记录
        L_record(end+1) = L_km;
        V_record(end+1) = v_kmh;
        H_record(end+1) = h;        
        % 找最大        
if L_km > max_L
            max_L = L_km;
            best_V = v_kmh;
            best_h = h;
            % 找到最佳点在数据表中的行号
            best_idx = find(data.feixinggaodu == h & data.xunhangsudu == v_kmh & data.dimianzaosheng <= noise_limit, 1);        
end    
end
end
if isnan(best_V)    
error('所有高度区间都没有满足噪声限制的速度组合！');
end

%------ 3. 计算总时间 ------

best_V_ms = best_V / 3.6;
t_cruise = (max_L * 1000) / best_V_ms;    % 巡航时间(s)
t_up = best_h / v_climb;                  % 起飞时间(s)
t_down = best_h / v_climb;                % 降落时间(s)
t_total = (t_cruise + t_up + t_down) / 3600; % 小时

%------ 4. 输出结果 -------

fprintf('------ 霄驰V1 满载最大航程（考虑地面噪声≤%.0fdB(A)）------\n', noise_limit);
fprintf('总质量：%.0f kg，最优巡航高度：%.0f m\n', m, best_h);
fprintf('最优巡航速度：%.2f km/h (%.2f m/s)\n', best_V, best_V_ms);
fprintf('最大航程：%.2f km\n', max_L);
fprintf('总飞行时间：%.2f 小时\n', t_total);
fprintf('最优点地面噪声：%.1f dB(A)\n', data.dimianzaosheng(best_idx));
fprintf('起飞+降落总能耗：%.2f kWh (%.2f%%电池)\n', ...    
(m * g * best_h / eta + 0.4 * m * g * best_h / eta)/3600/1000, ...    
(m * g * best_h / eta + 0.4 * m * g * best_h / eta)/E_bat_J*100);

%------ 5. 可视化 ------

figure;
scatter3(H_record, V_record, L_record, 40, L_record, 'filled');
xlabel('巡航高度 (m)');
ylabel('巡航速度 (km/h)');
zlabel('最大航程 (km)');
title('不同高度和速度下最大航程（满足噪声限制）');
colorbar;
grid on;

%------ 6. 能耗结构饼图 ------

E_up = m * g * best_h / eta;
E_down = 0.4 * E_up;
E_cruise = 0.5 * rho * S * Cd * best_V_ms^2 * (max_L*1000) / eta; % 巡航能耗(J)
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