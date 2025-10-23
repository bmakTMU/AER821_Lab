%% Using Collected Data from Lab

givendata = load('CMG_Part1_Wheel_speedkp0_14_Ki0_35_Kd1_6');

T = givendata.wheel_speed';

deg = unique(T(:,5));

A = zeros([size(T)-[0 1] size(deg)]);

for i = 1:length(deg)

    [n,m] = find(T==deg(i));
    for j =1:length(n)
        A(j,:,i) = T(n(j), 1:end-1);
    end
end
hold on
% Plot the data for each unique degree
for i = 1:length(deg)
    figure(i)

    avg = mean(A(:,1,i));
    stddev = std(A(:,1,i));
    xaxis = [max(A(:,1,i))-50,max(A(:,1,i))+10];

    subplot(3,1,1);
    plot(A(:,1,i),A(:,2,i));
    title(['Data for Degree: ', num2str(deg(i))]);
    xlabel('Time');
    ylabel('Degrees');
    xlim(xaxis);

    subplot(3,1,2);
    plot(A(:,1,i), A(:,3,i));
    xlabel('Time');
    ylabel('RPM');
    xlim(xaxis);

    subplot(3,1,3);
    plot(A(:,1,i),A(:,4,i));
    xlabel('Time');
    ylabel('Gyro Deflection Angle');
    xlim(xaxis);
    
    
    
end

hold off