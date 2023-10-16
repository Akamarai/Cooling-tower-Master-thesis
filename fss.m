function dwdTw = fss(w, iss, Tw, Ta, d)

    if d.a == false
        [psatw, wsa, cpa, cpv, cpw, ~] = checkss(w, Ta, d);
    else
        % Computation of pressure of water vapor evaluated at Tdb
        z = 10.79586 * (1 - 273.15 / Tw) + 5.02808 * log10(273.15 / Tw) ...
            +1.50474e-4 * (1 - (10 ^ (-8.29692 * (Tw / 273.15)) - 1)) ...
            +4.2873e-4 * (10 ^ (4.76955 * (1 - 273.15 / Tw)) - 1) ...
            + 2.786118312;

        psatw = 10 ^ z;

        % Computation of humidity ratio for saturated air
        wsa = ((2501.6 - 2.3263 * (Ta - 273.15)) / (2501.6 + 1.8577 * (Ta - 273.15) - 4.184 * (Ta - 273.15))) * ((0.62509 * psatw) / (d.P - 1.005 * psatw));

        [cpa, cpv, cpw] = heatcap(Tw);
    end

    % Calculate the humidity mass ratio at the gas-phase side at water temperature
    wsw = (0.6250 * psatw) / (d.P - 1.005 * psatw);

    % Enthalpy of water vapor at the bulk water temperature Tw [J/kg]
    iv = d.ifgwo + cpv * (Tw - 273.15);

    % Enthlpy of air-water vapor mixture [J/kg dry air K]
    imasw = cpa * (Tw - 273.15) + iv * wsw;

    %Water/Air ratio flow rate
    mwma = (d.mw / d.madry) * (1 - (d.madry / d.mw) * (d.wend - w));

    % Lewis factor
    Lef = 0.865 ^ (0.667) * ((((wsw + 0.622) / (wsa + 0.622)) - 1) / (log((wsw + 0.622) / (wsa + 0.622))));

    dwdTw = (cpw * mwma * (wsw - wsa)) / (imasw - iss + (Lef - 1) * ((imasw - iss - (wsw - wsa) * iv) + (w - wsa) * cpw * Tw) - (w - wsw) * cpw * Tw);
end
