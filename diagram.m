%% diagram
% harvesting model
r_h = 1;
K = 10;
harvest = @(X) r_h*X*(1-X/K)-c*X^2/(X^2+1);    

% Eutrophication model
a = 0.5;
h_e = 1;
eutrophication = @(X) a - r_e*X+c*X^8/(X^8+1);    

% Vegetation?turbidity model
h_v = 0.2;
r_v = 2;
veg_turb = @(X) 0.5*X*(1-X*(r_v^4+E^4)/r_v^4);    

%% harvesting bifurcation diagram
bifur_h = [];
c_range = 1:0.02:3;
for c = c_range
    syms X
    sol_h = vpasolve(r_h*X*(1-X/K)-c*X^2/(X^2+1), X);
    bifur_h = [bifur_h double(sol_h)];
end

figure;
hold on
plot(c_range, bifur_h(1,:))
plot(c_range(1:81), [bifur_h(2,1:40) bifur_h(4,41:81)])
plot(c_range(41:end), bifur_h(2,41:end))
plot(c_range(41:81), bifur_h(3,41:81))
hold off

%% Eutrophication bifurcation diagram
bifur_e = [];
e_range = 0.5:0.01:2;
for c = e_range
    syms X
    sol_e = vpasolve(a - h_e*X+c*X^8/(X^8+1), X);
    bifur_e = [bifur_e double(sol_e)];
end

figure;
hold on
plot(e_range, bifur_e(1,:))
plot(e_range(38:end), bifur_e(2,38:end))
plot(e_range(38:end), bifur_e(3,38:end))
hold off

%% Vegetation?turbidity model
bifur_v = [];
v_range = 2:0.1:12;
h_v = 0.2;
r_v = 2;
for c = v_range
    syms X
    sol_v = vpasolve(0.5*X*(1-X*(r_v^4+(c*h_v/(h_v+X))^4)/r_v^4), X);
    bifur_v = [bifur_v double(sol_v)];
end

figure;
hold on
plot(v_range, bifur_v(1,:))
plot(v_range(1:53), [bifur_v(2,1:32) bifur_v(4,33:53)])
plot(v_range(33:end), bifur_v(2,33:end))
plot(v_range(33:53), bifur_v(3,33:53))
hold off

%% plot
figure
subplot(1,3,1)
hold on
plot(c_range, bifur_h(1,:))
plot(c_range(1:81), [bifur_h(2,1:40) bifur_h(4,41:81)])
plot(c_range(41:end), bifur_h(2,41:end))
plot(c_range(41:81), bifur_h(3,41:81))
hold off
subplot(1,3,2)
hold on
plot(e_range, bifur_e(1,:))
plot(e_range(38:end), bifur_e(2,38:end))
plot(e_range(38:end), bifur_e(3,38:end))
hold off
subplot(1,3,3)
hold on
plot(v_range, bifur_v(1,:))
plot(v_range(1:53), [bifur_v(2,1:32) bifur_v(4,33:53)])
plot(v_range(33:end), bifur_v(2,33:end))
plot(v_range(33:53), bifur_v(3,33:53))
hold off