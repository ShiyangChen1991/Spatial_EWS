%% harvest simulation
load('r.mat');
size = 20;
sigma = 0.1;

[S_harvest, T] = run_simulation('harvest', size, [2, 1/16000], sigma, r, 600000, 0.02);
pop_harvest = [];
for i = 1:length(T)
     pop_harvest(i) = sum(S_harvest(i,:));
end
figure;
plot(2 + 1/16000*T,pop_harvest)

csvwrite('harvest.csv', S_harvest(100000:10:400000, :))
csvwrite('harvest_steady.csv', S_harvest(300000, :))

%% harvest comparison
harvest_steady = csvread('harvest_steady.csv');
harvest_shrink = csvread('harvest_shrink.csv', 1);
harvest_shrink_per = csvread('harvest_shrink_per.csv', 1);
S_harvest = csvread('harvest.csv' ,1);

harvest_steady = reshape(harvest_steady, [size, size]);
harvest_c_range = 2.1:0.01:2.45;
G = sigma^2*eye(size*size, size*size);
harvest_analytical = []; harvest_analytical_per = [];

for i = 1:length(harvest_c_range)
    [~, D_c] = analytical('harvest', harvest_c_range(i), r, harvest_steady, G);
    harvest_analytical(i) = max(diag(D_c));
    harvest_analytical_per(i) = max(diag(D_c))/norm(diag(D_c));
end

harvest_c_range_new = 2.1+1500/800*0.01:1000/800*0.01:2.1+1500/800*0.01+1000/800*0.01*26;

figure;
hold on
plot(harvest_c_range, harvest_analytical,'r--', 'linewidth',3)
plot(harvest_c_range_new(1:25), harvest_shrink(1:25), 'linewidth',2)
legend('analytical', 'estimation')
hold off

figure;
hold on
plot(harvest_c_range, harvest_analytical_per,'r--', 'linewidth',3)
plot(harvest_c_range_new(1:25), harvest_shrink_per(1:25), 'linewidth',2)
legend('analytical', 'estimation')
hold off
%% eutrophication simulation
load('r.mat');
size = 20;
sigma = 0.05;

[S_eutrophication, T] = run_simulation('eutrophication', size, [1.3, -1/16000], sigma, r, 600000, 0.02);
pop_eutrophication = [];
for i = 1:length(T)
     pop_eutrophication(i) = sum(S_eutrophication(i,:));
end
figure;
plot(1.3 - 1/16000*T,pop_eutrophication)

csvwrite('eutrophication.csv', S_eutrophication(100000:10:400000, :))
csvwrite('eutrophication_steady.csv', S_eutrophication(200000, :))

%% eutrophication comparison
eutrophication_steady = csvread('eutrophication_steady.csv');
eutrophication_shrink = csvread('eutrophication_shrink.csv', 1);
eutrophication_shrink_per = csvread('eutrophication_shrink_per.csv', 1);
S_eutrophication = csvread('eutrophication.csv' ,1);

eutrophication_steady = reshape(eutrophication_steady, [size, size]);
eutrophication_c_range = 0.85:0.02:1.2;
G = sigma^2*eye(size*size, size*size);
eutrophication_analytical = [];eutrophication_analytical_per = [];

for i = 1:length(eutrophication_c_range)
    [~, D_c] = analytical('eutrophication', eutrophication_c_range(i), r, eutrophication_steady, G);
    eutrophication_analytical(i) = max(diag(D_c));
    eutrophication_analytical_per(i) = max(diag(D_c))/norm(diag(D_c));
end

eutrophication_c_range_new = 1.2-1500/800*0.01:-1000/800*0.01:1.2-1500/800*0.01-1000/800*0.01*26;

figure;
hold on
plot(eutrophication_c_range, eutrophication_analytical,'r--', 'linewidth',3)
plot(eutrophication_c_range_new(1:24), eutrophication_shrink(1:24), 'linewidth',2)
legend('analytical', 'estimation')
hold off

figure;
hold on
plot(eutrophication_c_range, eutrophication_analytical_per,'r--', 'linewidth',3)
plot(eutrophication_c_range_new(1:24), eutrophication_shrink_per(1:24), 'linewidth',2)
legend('analytical', 'estimation')
hold off
%% Vegetation?turbidity simulation
load('r.mat');
size = 20;
sigma = 0.1;

[S_veg, T] = run_simulation('veg_turb', size, [2, 1/2000], sigma, 2*r, 600000, 0.02);
pop_veg_turb = [];
for i = 1:length(T)
     pop_veg_turb(i) = sum(S_veg(i,:));
end
figure;
plot(2 + 1/2000*T,pop_veg_turb)

csvwrite('veg_turb.csv', S_veg(100000:10:450000, :))
csvwrite('veg_turb_steady.csv', S_veg(300000, :))

%% Vegetation?turbidity comparison
veg_turb_steady = csvread('veg_turb_steady.csv');
veg_turb_shrink = csvread('veg_turb_shrink.csv', 1);
veg_turb_shrink_per = csvread('veg_turb_shrink_per.csv', 1);
S_veg = csvread('veg_turb.csv',1);

veg_turb_steady = reshape(veg_turb_steady, [size, size]);
veg_turb_c_range = 3:0.1:6.5;
G = sigma^2*eye(size*size, size*size);
veg_turb_analytical = [];

for i = 1:length(veg_turb_c_range)
    [~, D_c] = analytical('veg_turb', veg_turb_c_range(i), 2*r, veg_turb_steady, G);
    veg_turb_analytical(i) = max(diag(D_c));
    veg_turb_analytical_per(i) = max(diag(D_c))/norm(diag(D_c));
end

veg_turb_c_range_new = 3+1500/100*0.01:1000/100*0.01:3+1500/100*0.01+1000/100*0.01*31;

figure;
hold on
plot(veg_turb_c_range, veg_turb_analytical,'r--', 'linewidth',3)
plot(veg_turb_c_range_new(1:32), veg_turb_shrink(1:32), 'linewidth',2)
legend('analytical', 'estimation')
hold off

figure;
hold on
plot(veg_turb_c_range, veg_turb_analytical_per,'r--', 'linewidth',3)
plot(veg_turb_c_range_new(1:32), veg_turb_shrink_per(1:32), 'linewidth',2)
legend('analytical', 'estimation')
hold off
%% plot
figure
subplot(3,3,1)
plot(2 + 1/16000*T,pop_harvest)
subplot(3,3,2)
plot(1.3 - 1/16000*T,pop_eutrophication)
subplot(3,3,3)
plot(2 + 1/2000*T,pop_veg_turb)
subplot(3,3,4)
hold on
plot(harvest_c_range, harvest_analytical,'r--', 'linewidth',3)
plot(harvest_c_range_new(1:24), harvest_shrink(1:24), 'linewidth',2)
legend('analytical', 'estimation')
hold off
subplot(3,3,5)
hold on
plot(eutrophication_c_range, eutrophication_analytical,'r--', 'linewidth',3)
plot(eutrophication_c_range_new(1:23), eutrophication_shrink(1:23), 'linewidth',2)
legend('analytical', 'estimation')
hold off
subplot(3,3,6)
hold on
plot(veg_turb_c_range, veg_turb_analytical,'r--', 'linewidth',3)
plot(veg_turb_c_range_new(1:32), veg_turb_shrink(1:32), 'linewidth',2)
legend('analytical', 'estimation')
hold off
subplot(3,3,7)
hold on
plot(harvest_c_range, harvest_analytical_per,'r--', 'linewidth',3)
plot(harvest_c_range_new(1:24), harvest_shrink_per(1:24), 'linewidth',2)
legend('analytical', 'estimation')
hold off
subplot(3,3,8)
hold on
plot(eutrophication_c_range, eutrophication_analytical_per,'r--', 'linewidth',3)
plot(eutrophication_c_range_new(1:23), eutrophication_shrink_per(1:23), 'linewidth',2)
legend('analytical', 'estimation')
hold off
subplot(3,3,9)
hold on
plot(veg_turb_c_range, veg_turb_analytical_per,'r--', 'linewidth',3)
plot(veg_turb_c_range_new(1:32), veg_turb_shrink_per(1:32), 'linewidth',2)
legend('analytical', 'estimation')
hold off