clc 
clear all

%% Load data 

% 5e-3 Nm
load current5.mat 
load current_torque5.mat
load current_volts5.mat

load speed5.mat
load speed_torque5.mat
load speed_voltage5.mat

current5 = current';
currentVoltage5 = currentVoltage';
currentTorqu5 = currentTorqu';

speed5 = speed';
speedVoltage5 = speedVoltage';
speedTorque5 = speedTorque';

% 10e-3 Nm
load current10.mat 
load current_torque10.mat
load current_volts10.mat

load speed10.mat
load speed_torque10.mat
load speed_voltage10.mat

current10 = current';
currentVoltage10 = currentVoltage';
currentTorqu10 = currentTorqu';

speed10 = speed';
speedVoltage10 = speedVoltage';
speedTorque10 = speedTorque';

%% Current plots

% Current
figure(1) 
subplot(2,1,1)
hold on
for i=2:width(current5)
    plot(current5(:,1), current5(:,i), 'LineWidth', 1.5);
end 
xlabel('Time');
ylabel('Current (A)');
title('Time vs Current, 5 mNm');
legend('Commanded Torque', "Load Cell Torque", "Computed Torque (from current)")
grid on; hold off;


subplot(2,1,2)
hold on
for i=2:width(current10)
    plot(current10(:,1), current10(:,i), 'LineWidth', 1.5);
end 
xlabel('Time');
ylabel('Current (A)');
title('Time vs Current, 10 mNm');
legend('Commanded Torque', "Load Cell Torque", "Computed Torque (from current)")
grid on; hold off;


% Voltage
figure(2)
subplot(2,1,1)
hold on
plot(currentVoltage5(:,1), currentVoltage5(:,2), 'LineWidth', 1.5);
xlabel('Time');
ylabel('Voltage (V)');
title('Time vs Voltage, 5 mNm');
grid on; hold off;

subplot(2,1,2)
hold on
plot(currentVoltage10(:,1), currentVoltage10(:,2), 'LineWidth', 1.5);
xlabel('Time');
ylabel('Voltage (V)');
title('Time vs Voltage, 10 mNm');
grid on; hold off;

% Torque

figure(3) 
subplot(2,1,1)
hold on
for i=2:width(currentTorqu5)
    plot(currentTorqu5(:,1), currentTorqu5(:,i), 'LineWidth', 1.5);
end 
xlabel('Time');
ylabel('Current (A)');
title('Time vs Current, 5 mNm');
legend('Commanded Torque', "Load Cell Torque", "Computed Torque (from current)")
grid on; hold off;


subplot(2,1,2)
hold on
for i=2:width(currentTorqu10)
    plot(currentTorqu10(:,1), currentTorqu10(:,i), 'LineWidth', 1.5);
end 
xlabel('Time');
ylabel('Current (A)');
title('Time vs Current, 10 mNm');
legend('Commanded Torque', "Load Cell Torque", "Computed Torque (from current)")
grid on; hold off;


% 
% figure(3)
% subplot(2,1,1)
% hold on
% plot(currentTorqu5(:,1), currentTorqu5(:,2), 'LineWidth', 1.5);
% xlabel('Time');
% ylabel('Torque (Nm)');
% title('Time vs Torque, 5 mNm');
% grid on; hold off;
% 
% subplot(2,1,2)
% hold on
% plot(currentTorqu10(:,1), currentTorqu10(:,2), 'LineWidth', 1.5);
% xlabel('Time');
% ylabel('Torque (Nm)');
% title('Time vs Torque, 10 mNm');
% grid on; hold off;




%% Speed plots

% Current
figure(4) 
subplot(2,1,1)
hold on
for i=2:width(speed5)
    plot(speed5(:,1), speed5(:,i), 'LineWidth', 1.5);
end 
xlabel('Time');
ylabel('Current (A)');
title('Time vs Current, 5 mNm');
legend('Commanded Torque', "Load Cell Torque", "Computed Torque (from speed)")
grid on; hold off;


subplot(2,1,2)
hold on
for i=2:width(speed10)
    plot(speed10(:,1), speed10(:,i), 'LineWidth', 1.5);
end 
xlabel('Time');
ylabel('Current (A)');
title('Time vs Current, 10 mNm');
legend('Commanded Torque', "Load Cell Torque", "Computed Torque (from current)")
grid on; hold off;


% Voltage
figure(5)
subplot(2,1,1)
hold on
plot(speedVoltage5(:,1), speedVoltage5(:,2), 'LineWidth', 1.5);
xlabel('Time');
ylabel('Voltage (V)');
title('Time vs Voltage, 5 mNm');
grid on; hold off;

subplot(2,1,2)
hold on
plot(speedVoltage10(:,1), speedVoltage10(:,2), 'LineWidth', 1.5);
xlabel('Time');
ylabel('Voltage (V)');
title('Time vs Voltage, 10 mNm');
grid on; hold off;

% Torque
figure(6)
subplot(2,1,1)
hold on
plot(speedTorque5(:,1), speedTorque5(:,2), 'LineWidth', 1.5);
xlabel('Time');
ylabel('Torque (Nm)');
title('Time vs Torque, 5 mNm');
grid on; hold off;

subplot(2,1,2)
hold on
plot(speedTorque10(:,1), speedTorque10(:,2), 'LineWidth', 1.5);
xlabel('Time');
ylabel('Torque (Nm)');
title('Time vs Torque, 10 mNm');
grid on; hold off;
