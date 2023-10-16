function [cpa, cpv, cpw] = heatcap(Temp)

    T = (Temp + 273.15) / 2;
    cpa = 1.045356e3 -3.161783e-1 * T +7.083814e-4 * T ^ 2 -2.705209e-7 * T ^ 3;
    % Specific heat of saturated water wapor [J/kgK]
    cpv = 1.3605e3 + 2.31334 * T -2.46784e-10 * T ^ 5 +5.91332e-13 * T ^ 6;
    % Specific heat of water [J/kgK]
    cpw = 8.15599e3 - 2.80627 * 10 * T +5.11283e-2 * T ^ 2 -2.17582e-13 * T ^ 6;

end
