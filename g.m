function dimadTw = g(w, ima, Tw, d)

    [cpa, cpv, cpw] = heatcap(Tw);

    % Computation of pressure of water vapor evaluated at Tw
    z = 10.79586 * (1 - 273.15 / Tw) + 5.02808 * log10(273.15 / Tw) ...
        +1.50474e-4 * (1 - (10 ^ (-8.29692 * (Tw / 273.15)) - 1)) ...
        +4.2873e-4 * (10 ^ (4.76955 * (1 - 273.15 / Tw)) - 1) ...
        + 2.786118312;

    psatw = 10 ^ z;

    % Calculate the humidity mass ratio at the gas-phase side at water temperature
    wsw = (0.6250 * psatw) / (d.P - 1.005 * psatw);

    % Enthalpy of water vapor at the bulk water temperature Tw [J/kg]
    iv = d.ifgwo + cpv * (Tw - 273.15);

    % Enthlpy of air-water vapor mixture [J/kg dry air K]
    imasw = cpa * (Tw - 273.15) + iv * wsw;

    %Water/Air ratio flow rate
    mwma = (d.mw / d.madry) * (1 - (d.madry / d.mw) * (d.wend - w));

    % Lewis factor
    Lef = 0.865 ^ (0.667) * ((((wsw + 0.622) / (w + 0.622)) - 1) / (log((wsw + 0.622) / (w + 0.622))));

    dimadTw = mwma * cpw * (1 + ((wsw - w) * cpw * Tw) / (imasw - ima + (Lef - 1) * (imasw - ima - (wsw - w) * iv) - (wsw - w) * cpw * Tw));

end
