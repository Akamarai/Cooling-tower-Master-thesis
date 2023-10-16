function [pv_wb, wsa, cpa, cpv, cpw, iss] = checkss(w, Ta, d)

    % Computation of pressure of water vapor evaluated at Tdb
    z = 10.79586 * (1 - 273.15 / Ta) + 5.02808 * log10(273.15 / Ta) ...
        +1.50474e-4 * (1 - (10 ^ (-8.29692 * (Ta / 273.15)) - 1)) ...
        +4.2873e-4 * (10 ^ (4.76955 * (1 - 273.15 / Ta)) - 1) ...
        + 2.786118312;

    pv_wb = 10 ^ z;

    % Computation of humidity ratio for saturated air
    wsa = ((2501.6 - 2.3263 * (Ta - 273.15)) / (2501.6 + 1.8577 * (Ta - 273.15) - 4.184 * (Ta - 273.15))) * ((0.62509 * pv_wb) / (d.P - 1.005 * pv_wb));

    [cpa, cpv, cpw] = heatcap(Ta);

    iss = cpa * (Ta - 273.15) + wsa * (d.ifgwo + cpv * (Ta - 273.15)) + (w - wsa) * cpw * (Ta - 273.15);
end
