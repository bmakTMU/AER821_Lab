 givendata = load('CMG_Part1_Wheel_speedkp0_14_Ki0_35_Kd1_6');
T = givendata.wheel_speed';

deg = unique(T(:,5));

% Preallocate cell array (easier than forcing into 3D matrix)
A = cell(length(deg),1);

for i = 1:length(deg)
    % Only select rows where column 5 (command angle) equals current degree
    rows = T(:,5) == deg(i);
    A{i} = T(rows, 1:end-1);   % store all but the degree column
end

hold on
for i = 1:length(deg)
    figure(i)
    clf

    data = A{i};
    time = data(:,1);
    actual = data(:,2);      % actual response (angle)
    rpm = data(:,3);
    gyro = data(:,4);
    desired = deg(i) * ones(size(time));

    xaxis = [min(time), max(time)];

    % ---- Angle (Actual vs Desired) ----
    subplot(3,1,1);
    plot(time, actual, 'b', 'LineWidth', 1.3); hold on;
    plot(time, desired, 'r--', 'LineWidth', 1.4);
    title(['Response for Command Angle = ', num2str(deg(i)), '°']);
    xlabel('Time (s)');
    ylabel('Angle (°)');
    legend('Actual Output', 'Desired Output');
    grid on; xlim(xaxis);
end
hold off