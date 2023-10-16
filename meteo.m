% Read the data from the Excel file, including date-time, temperature, and humidity
data = readcell("meteo.xlsx");

% Calculate the average temperature for each month
months = month(data(: , 1));
uniqueMonths = unique(months);
monthlyAverageTemperatures = zeros(size(uniqueMonths));

for i = 1:length(uniqueMonths)
    monthIndex = uniqueMonths(i);
    monthlyAverageTemperatures(i) = mean(temperatureData(months == monthIndex));
end

% Create a bar plot to display the monthly average temperatures
bar(uniqueMonths, monthlyAverageTemperatures);
xlabel('Month');
ylabel('Average Temperature');
title('Monthly Average Temperature for 2023');
xticks(1:12);
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'});

% Optionally save the plot as an image file
% saveas(gcf, 'monthly_average_temperature.png');
