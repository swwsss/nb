% 杭电数模赛A题 问题四主流程
% 1. 读入九站经纬度，计算站点间距离矩阵
% 2. 给定任意一条路径，自动计算每段能耗/所需最小电量/充电时间
% 3. 模拟全程飞行+充电，输出各段参数，推荐两条方案

clc; clear;

%% 1. 站点及经纬度
station_names = {'萧山机场','钱塘金沙湖','西湖文化广场','上城吴山广场',...
    '临安人民广场','富阳秦望广场','桐庐中心广场','建德新安江广场','淳安千岛湖广场'};
locations = [
    30.2295, 120.1699;   % 萧山机场站
    30.3111, 120.3458;   % 钱塘金沙湖站
    30.2765, 120.1654;   % 西湖文化广场站
    30.2412, 120.1705;   % 上城吴山广场站
    30.2333, 119.7147;   % 临安人民广场站
    30.0486, 119.9508;   % 富阳秦望广场站
    29.7975, 119.6736;   % 桐庐中心广场站
    29.4758, 119.2803;   % 建德新安江广场站
    29.6089, 119.0225;   % 淳安千岛湖广场站
];
num_sites = size(locations, 1);

%% 2. 距离矩阵
distance_matrix = zeros(num_sites);
for i = 1:num_sites
    for j = 1:num_sites
        if i ~= j
            [distance_deg, ~] = distance(locations(i,1), locations(i,2), ...
                                         locations(j,1), locations(j,2), ...
                                         'degrees');
            distance_m = deg2km(distance_deg) * 1000; % m
            distance_matrix(i,j) = distance_m;
        end
    end
end
distance_km = round(distance_matrix / 1000, 2);

disp('站点间距离矩阵（单位：公里）:');
disp(array2table(distance_km, 'VariableNames', station_names, 'RowNames', station_names));

%% 3. eVTOL参数（可根据前几问定值）
g = 9.79; rho = 1.18; S = 2; Cd = 0.3; eta = 0.8;
bat_mass = 200; bat_den = 300;
E_bat_Wh = bat_mass * bat_den;
E_bat_J = E_bat_Wh * 3600;
total_mass = 800 + 5 * 65 + 4*20; % 满载1300kg
cruise_v = 140 / 3.6;  % 140km/h为例，可调整
h_cruise = 400;        % 400m为例，可调整
v_climb = 5;           % m/s

% 充电参数
charge_rate = 0.1 * E_bat_J / 60;  % 10%/min（Joule/min）
min_charge_time = 1; % 最少充电1分钟

%% 4. 路线规划（可手动设定两条方案，至少5站，含机场、千岛湖、吴山广场）
% 推荐线路1：最短距离常规线路
route1 = [1,2,4,6,7,9]; % 萧山机场→金沙湖→吴山广场→富阳→桐庐→千岛湖
% 推荐线路2：中间多停靠、充电次数较少
route2 = [1,3,4,5,8,9]; % 萧山机场→西湖文化广场→吴山广场→临安→建德→千岛湖

routes = {route1, route2};
num_routes = length(routes);

%% 5. 主流程：依次模拟每条线路
for route_id = 1:num_routes
    route = routes{route_id};
    fprintf('\n---- 推荐线路%d ----\n', route_id);
    fprintf('站点顺序: ');
    disp(strjoin(station_names(route), ' → '));
    
    energy_left = E_bat_J; % 每条线路起点电量100%
    energy_log = zeros(1, length(route));
    charge_time_log = zeros(1, length(route));
    total_flight_time = 0;
    total_charge_time = 0;
    
    for seg = 1:length(route)-1
        i = route(seg);
        j = route(seg+1);
        dist = distance_matrix(i,j);
        
        % 能耗计算
        E_up = total_mass * g * h_cruise / eta;
        E_down = 0.4 * E_up;
        E_cruise = 0.5 * rho * S * Cd * cruise_v^2 * dist / eta;
        E_flight = E_up + E_down + E_cruise;
        t_cruise = dist / cruise_v;
        t_up = h_cruise / v_climb;
        t_down = h_cruise / v_climb;
        t_flight = t_cruise + t_up + t_down;
        
        % 如果是上城吴山广场，到站后要观光2小时，消耗30%电量
        if j == 4
            E_sight = 0.3 * E_bat_J;
            t_sight = 2*3600; % 2小时
            fprintf('到%s站观光2小时，额外消耗电量%.2f%%。\n', station_names{j}, 30);
        else
            E_sight = 0;
            t_sight = 0;
        end
        
        % 判断剩余电量是否足够本段飞行+观光+30%冗余
        E_need = E_flight + E_sight + 0.3 * E_bat_J;
        if energy_left < E_need
            % 需充电
            E_charge = min(E_bat_J, E_need) - energy_left;
            charge_time = max(min_charge_time, ceil(E_charge / charge_rate)); % 分钟
            energy_left = energy_left + charge_time * charge_rate;
            if energy_left > E_bat_J
                energy_left = E_bat_J; % 不超过满电
            end
            total_charge_time = total_charge_time + charge_time;
            charge_time_log(seg) = charge_time;
            fprintf('在%s站充电%d分钟（补%.1f%%电量）\n', station_names{i}, charge_time, 100*E_charge/E_bat_J);
        else
            charge_time_log(seg) = 0;
        end
        
        % 执行飞行+观光
        energy_left = energy_left - (E_flight + E_sight);
        total_flight_time = total_flight_time + (t_flight + t_sight)/60; % 分钟
        if energy_left < 0
            fprintf('能量不足，航段失败！\n');
            break
        end
        energy_log(seg+1) = energy_left;
        fprintf('%s→%s 航段%.2fkm，飞行耗时%.2f分钟，剩余%.2f%%电量\n', ...
            station_names{i}, station_names{j}, dist/1000, (t_flight+t_sight)/60, 100*energy_left/E_bat_J);
    end
    fprintf('本线路总飞行时间：%.1f分钟，总充电时间：%.1f分钟，共充电%d次。\n', ...
        total_flight_time, total_charge_time, sum(charge_time_log>0));
    fprintf('到达终点剩余电量：%.1f%%\n', 100*energy_left/E_bat_J);

    % 可视化剩余电量变化
    figure;
    plot(1:length(route), energy_log/E_bat_J*100, '-o','LineWidth',2);
    xticks(1:length(route));
    xticklabels(station_names(route));
    ylabel('剩余电量(%)');
    title(sprintf('推荐线路%d电量变化曲线', route_id));
    grid on;
end
