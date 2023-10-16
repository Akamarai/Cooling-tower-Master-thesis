close all
clear
clc

%% Input data test
data = readcell("tabella_dati.xlsx");

d = struct(); % Imput data della torre da validare

d.Tairin = round(cell2mat(data(2:end,3))); % °C
d.Tairindry = round(cell2mat(data(2:end,2))); % °C
d.Pair = cell2mat(data(2:end,5))*100; % Pa
d.Twout = 27.77; % °C
d.Twin = 50; % °C
d.mw = 12500; % kg/s
N = 5; % Number of steps of RK method

%% Loop
for f=1:length(d.Tairin)

    % Model function calling
    [Tdb, Twb, Tw, ima, w, Mep,water,d] = modellosupersaturo(d.Tairindry(f), d.Tairin(f),d.Pair(f), N, d);
    pTdb(f,:) = Tdb(1:end-1);
    pTwb(f,:) = Twb;
    pTw(f,:) = Tw;
    pima(f,:) = ima;
    pw(f,:) = w(1:end-1);
    pMep(f,:) = Mep;
    pwater(f,:) = water;
    % Create a table with the data
    dataTable=array2table([(Tdb(1:end-1)' -273.15), (Twb'-273.15), (Tw' -273.15), ima(1:end-1)', w(1:end-1)', Mep'],'VariableNames',{'Tdb','Twb','Tw','ima','w','Mep'});
    tabella = sprintf('data_iteration_%d.xlsx', f);
    writetable(dataTable, tabella);

    %% Computation of parameter

    % Convergence check
    ERRima(f)= abs((ima(end) - ima(end -1)) / (ima(end) + ima(end -1)));

    %Water flow rate evaporated
    mwevap(f) = d.madry * (w(end - 1) - w(1));
    mwevap_ct(f) = d.madry * (w(end) - w(1));

    %Heat calculation
    sensible_heat_transfer(f) = d.madry * (ima(end-1) - ima(1));
    cooling_power_evap(f) = mwevap(f) * d.ifgwo;
    cooling_power_total(f) = cooling_power_evap(f) + sensible_heat_transfer(f) ;

    if Twb(1) >= Tdb(1)
        DT1 = Tw(1) - Tdb(1);
        DT2 = Tw(end) - Tdb(end -1);
    else
        DT1 = Tdb(1) - Tw(1);
        DT2 = Tdb(end) - Tw(end -1);
    end

    LMTD(f) = abs((DT1 - DT2) / log(DT1/DT2));
    Q(f) = LMTD(f) * 4.186 * water(end -1)
    %% Create the plot

    % Optimization for Tdb
    x = 0:length(Tdb) - 2;
    coefficients = polyfit(x, Tdb(1:N+1), 2);
    x_interp = linspace(min(x), max(x), 100);
    Tdb_interp = polyval(coefficients, x_interp);

    % Optimization for Twb
    x4 = 0:length(Twb) -1;
    coefficients4 = polyfit(x4, Twb, 2);
    x4_interp = linspace(min(x), max(x), 100);
    Twb_interp = polyval(coefficients4, x4_interp);

    % Optimization for ima
    x2 = 0:length(ima) -2;
    coefficients2 = polyfit(x2, ima(1:N+1), 2);
    x2_interp = linspace(min(x2), max(x2), 100);
    ima_interp = polyval(coefficients2, x2_interp);

    % Optimization for w

    x3 = 0:length(w)-2;
    coefficients3 = polyfit(x3, w(1:N+1), 2);
    x3_interp = linspace(min(x3), max(x3), 100);
    w_interp = polyval(coefficients3, x3_interp);

    % Optimization for Mep
    x5 = 0:length(Mep) -1;
    coefficients5 = polyfit(x5, Mep, 2);
    x5_interp = linspace(min(x), max(x), 100);
    Mep_interp = polyval(coefficients5, x5_interp);


    figure(f);
    hold on;
    plot(0:N, Tw-273.15, 'LineWidth', 2); % Plot Tw
    plot(x_interp, Tdb_interp-273.15, 'LineWidth', 2); % Plot Tdb
    plot(x4_interp, Twb_interp-273.15, 'LineWidth', 2);
    hold off;
    xlabel('Level of Spray zone');
    ylabel('Temperature [K]');
    legend('Tw', 'Tdb', 'Twb');
    title('Variation of Tw and Tdb over iterations');
    filename = sprintf('Tw_Tdb_Twb_%d.png', f);
    saveas(gcf, filename);
    close(gcf);

    figure(f+1);
    hold on;
    plot(0:N, w(1:N+1), 'LineWidth', 2); % Plot w
    hold off;
    xlabel('Level of Spray zone');
    ylabel('[kg/kg dry]');
    title('Humidity mass ratio');
    filename = sprintf('Humidity_mass_ratio_%d.png', f);
    saveas(gcf, filename);
    close(gcf);

    figure(f+2);
    hold on;
    plot(0:N, ima(1:N+1), 'LineWidth', 2); % Plot ima
    plot(x2_interp, ima_interp, 'LineWidth', 2);
    hold off;
    xlabel('Level of Spray zone');
    ylabel('Enthalpy J/kg dry air');
    legend('ima');
    filename = sprintf('Enthalpy_of_Spray_zone_%d.png', f);
    saveas(gcf, filename);
    close(gcf);

    figure(f+3);
    hold on;
    plot(0:N, water, 'LineWidth', 2); % Plot mw
    hold off;
    xlabel('Level of Spray zone');
    ylabel('[kg/s]');
    legend('mw');
    title('Water flow rate');
    filename = sprintf('Water_flow_rate_%d.png', f);
    saveas(gcf, filename);
    close(gcf);

end

axhandle=psychplotting(1,40,1,40);
figure(1)
hold on
for a=1:f
    plot(axhandle,pTdb(a,:)-273.15,pw(a,:)*1000,'LineWidth',1.5);
end
hold off
saveas(gcf, 'PsychroChart.png');

figure(2)
hold on
for b=1:f
    plot(0:N,pTw(b,:)-273.15);
end
hold off
saveas(gcf, 'Watercomparison.png');
