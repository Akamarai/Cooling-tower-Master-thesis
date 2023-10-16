function [Tdb, Twb, Tw, ima, w, Mep, water, d] = modellosupersaturo(Tairindry, Tairin,Pair, N, d)
%% Inizialization of vectors
water=d.mw;
Tdb = zeros(1, N + 2);
Twb = zeros(1, N + 1);
Tw = zeros(1, N + 1);
w = zeros(1, N + 2);
ima = zeros(1, N + 2);
Mep = zeros(1, N + 1);
d.P = Pair;
%% Computation of Air properties temperature
% Tdb = Dry bulb K
% Twb = Wet buld K
% w   = Humidity ratio kg/kg
% ima = Enthalpy of gas phase J/kg dry air

Tdb(1) = Tairindry +273.15;
Twb(1) = Tairin +273.15;

% Computation of pressure of water vapor evaluated at Twb
z(1) = 10.79586 * (1 - 273.15 / Twb(1)) + 5.02808 * log10(273.15 / Twb(1)) ...
    +1.50474e-4 * (1 - (10 ^ (-8.29692 * (Twb(1) / 273.15)) - 1)) ...
    +4.2873e-4 * (10 ^ (4.76955 * (1 - 273.15 / Twb(1))) - 1) ...
    + 2.786118312;

pv_wb = 10 ^ z(1);

%Merkel Number
Mep(1) = 1.55;

% Alternative to calculate pv
% pv_wb_alt = IAPWS_IF97('psat_T', Tdb(1)) * 10^6;

% Computation of humidity ratio for saturated air
w(1) = ((2501.6 - 2.3263 * (Twb(1) - 273.15)) / (2501.6 + 1.8577 * (Twb(1) - 273.15) - 4.184 * (Twb(1) - 273.15))) * ...
    ((0.62509 * pv_wb) / (d.P - 1.005 * pv_wb)) - ...
    ((1.00416 * ((Tdb(1) - 273.15) - (Twb(1) - 273.15))) / (2501.6 + 1.8577 * (Tdb(1) - 273.15) - 4.184 * (Twb(1) - 273.15)));

d.wo = w(1); % Value to put in database for computation of mass flow rate in f,g,h

%% Computation of enthalpy and specific heat capacity

T = (Tdb(1) + 273.15) / 2; % Temperature used as reference for the computation

% Specific heat of dry air [J/kgK]
cpa = 1.045356e3 -3.161783e-1 * T +7.083814e-4 * T ^ 2 -2.705209e-7 * T ^ 3;

% Specific heat of saturated water wapor [J/kgK]
cpv = 1.3605e3 + 2.31334 * T -2.46784e-10 * T ^ 5 +5.91332e-13 * T ^ 6;

%Latent heat of water evaluated at 0°C
d.ifgwo = 3.4831814e6 -5.8627703e3 * 273.15 + 12.139568 * 273.15 ^ 2 -1.40290431e-2 * 273.15 ^ 3;

% Enthlpy of air-water vapor mixture [J/kg dry air K]
ima(1) = cpa * (Tdb(1) - 273.15) + w(1) * (d.ifgwo + cpv * (Tdb(1) - 273.15));
%% Initial approximation of variables
Tw(1) = ((d.Twin + 273.15) + 2 * Twb(1) + Tdb(1)) / 4;  % Outlet water temperature approximation based on Poppe Models [K]
Tdb(end) = ((d.Twin + 273.15 + Tw(1)) / 2); % Outlet air temperature approximation based on Poppe Models [K]
TEND = (d.Twin + 273.15 + Tw(1)) / 2; % Temperature used as reference for the end value of computation

%Specific heat of water [J/kgK] evaluated a (Twin+Twout)/2
cpw(1) = 8.15599e3 - 2.80627 * 10 * TEND +5.11283e-2 * TEND ^ 2 -2.17582e-13 * TEND ^ 6;

z(2) = 10.79586 * (1 - 273.15 / Tdb(end)) + 5.02808 * log10(273.15 / Tdb(end)) ...
    +1.50474e-4 * (10 ^ (-8.29692 * (Tdb(end) / 273.15) - 1)) ...
    +4.2873e-4 * (10 ^ (4.76955 * (1 - 273.15 / Tdb(end))) - 1) ...
    + 2.786118312;

pv_wb(2) = 10 ^ z(2);
w(end) = (0.6250 * pv_wb(2)) / (d.P - 1.005 * pv_wb(2));
d.wend = w(end);
T2 = (Tdb(end) + 273.15) / 2;

cpa(2) = 1.045356e3 -3.161783e-1 * T2 +7.083814e-4 * T2 ^ 2 -2.705209e-7 * T2 ^ 3;
cpv(2) = 1.3605e3 + 2.31334 * T2 -2.46784e-10 * T2 ^ 5 +5.91332e-13 * T2 ^ 6;
ima(end) = cpa(2) * (Tdb(end) - 273.15) + w(end) * (d.ifgwo + cpv(2) * (Tdb(end) - 273.15));
% Vapor mass flow rate [kg/s]
mair = (d.mw * cpw(1) * ((d.Twin + 273.15) - Tw(1))) / (ima(end) - ima(1));

% Dry air mass flow rate [kg dry/s]
d.madry = (2 * mair) / (2 + d.wo + d.wend);

%% Application of Runge-Kutta method to determine the working status of the cooling tower

DTW = ((d.Twin + 273.15) - Tw(1)) / N;
DELTA = 0.005;
d.a = false;

% Runge-Kutta 4th order to resolve the ODEs system of Cooling-Tower
for n = 1:N
    if d.a==false
        %normal cycle
        j(n + 1, 1) = DTW * f(w(n), ima(n), Tw(n), d);
        k(n + 1, 1) = DTW * g(w(n), ima(n), Tw(n), d);
        l(n + 1, 1) = DTW * h(w(n), ima(n), Tw(n), d);

        j(n + 1, 2) = DTW * f(w(n) + (j(n + 1, 1)) / 2, ima(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, d);
        k(n + 1, 2) = DTW * g(w(n) + (j(n + 1, 1)) / 2, ima(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, d);
        l(n + 1, 2) = DTW * h(w(n) + (j(n + 1, 1)) / 2, ima(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, d);

        j(n + 1, 3) = DTW * f(w(n) + (j(n + 1, 2)) / 2, ima(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, d);
        k(n + 1, 3) = DTW * g(w(n) + (j(n + 1, 2)) / 2, ima(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, d);
        l(n + 1, 3) = DTW * h(w(n) + (j(n + 1, 2)) / 2, ima(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, d);

        j(n + 1, 4) = DTW * f(w(n) + j(n + 1, 3), ima(n) + k(n + 1, 3), Tw(n) + DTW, d);
        k(n + 1, 4) = DTW * g(w(n) + j(n + 1, 3), ima(n) + k(n + 1, 3), Tw(n) + DTW, d);
        l(n + 1, 4) = DTW * h(w(n) + j(n + 1, 3), ima(n) + k(n + 1, 3), Tw(n) + DTW, d);

        w(n + 1) = w(n) + (j(n + 1, 1) + 2 * j(n + 1, 2) + 2 * j(n + 1, 3) + j(n + 1, 4)) / 6;
        ima(n + 1) = ima(n) + (k(n + 1, 1) + 2 * k(n + 1, 2) + 2 * k(n + 1, 3) + k(n + 1, 4)) / 6;
        Mep(n + 1) = Mep(n) + (l(n + 1, 1) + 2 * l(n + 1, 2) + 2 * l(n + 1, 3) + l(n + 1, 4)) / 6;

        % Calculate air temperature [K]
        Tw(n + 1) = Tw(n) + DTW;
        Tdb(n+1) = dry(w(n+1), Tw(n+1), ima(n+1), d);
        [~, ~, ~, ~, ~, ~, wet] = Psychrometricsnew('Tdb',(Tdb(n+1)-273.15),'w',w(n+1));
        Twb(n+1) = wet + 273.15;
        [~, ~, ~, ~, ~, iss(n + 1)] = checkss(w(n + 1), Tdb(n + 1), d);
        imact = checkima(w(n + 1), Tdb(n + 1), d);
        ERRima = abs((imact - iss(n + 1)) / (imact + iss(n + 1)));
    end

    if n > 1 && Twb(n) > Tdb(n) && d.a == false && ERRima > DELTA
        [~, ~, ~, ~, ~, iss(n)] = checkss(w(n), Tdb(n), d);
        X = ['Aria sovrassatura, ciclo n= ', num2str(n)];
        disp(X);
        d.a = true;
        %Supersaturated equations
        j(n + 1, 1) = DTW * fss(w(n), iss(n), Tw(n), Tdb(n), d);
        k(n + 1, 1) = DTW * gss(w(n), iss(n), Tw(n), Tdb(n), d);
        l(n + 1, 1) = DTW * hss(w(n), iss(n), Tw(n), Tdb(n), d);

        j(n + 1, 2) = DTW * fss(w(n) + (j(n + 1, 1)) / 2, iss(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
        k(n + 1, 2) = DTW * gss(w(n) + (j(n + 1, 1)) / 2, iss(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
        l(n + 1, 2) = DTW * hss(w(n) + (j(n + 1, 1)) / 2, iss(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, Tdb(n), d);

        j(n + 1, 3) = DTW * fss(w(n) + (j(n + 1, 2)) / 2, iss(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
        k(n + 1, 3) = DTW * gss(w(n) + (j(n + 1, 2)) / 2, iss(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
        l(n + 1, 3) = DTW * hss(w(n) + (j(n + 1, 2)) / 2, iss(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, Tdb(n), d);

        j(n + 1, 4) = DTW * fss(w(n) + j(n + 1, 3), iss(n) + k(n + 1, 3), Tw(n) + DTW, Tdb(n), d);
        k(n + 1, 4) = DTW * gss(w(n) + j(n + 1, 3), iss(n) + k(n + 1, 3), Tw(n) + DTW, Tdb(n), d);
        l(n + 1, 4) = DTW * hss(w(n) + j(n + 1, 3), iss(n) + k(n + 1, 3), Tw(n) + DTW, Tdb(n), d);

        w(n + 1) = w(n) + (j(n + 1, 1) + 2 * j(n + 1, 2) + 2 * j(n + 1, 3) + j(n + 1, 4)) / 6;
        ima(n + 1) = iss(n) + (k(n + 1, 1) + 2 * k(n + 1, 2) + 2 * k(n + 1, 3) + k(n + 1, 4)) / 6;
        Mep(n + 1) = Mep(n) + (l(n + 1, 1) + 2 * l(n + 1, 2) + 2 * l(n + 1, 3) + l(n + 1, 4)) / 6;

        % Calculate air temperature [K]
        Tw(n + 1) = Tw(n) + DTW;
        % Tdb(n+1) = dry(w(n+1), Tw(n+1), ima(n+1), d);
        [~, ~, ~, ~, ~, ~, wet] = Psychrometricsnew('h',ima(n+1),'w',w(n+1));
        Twb(n+1) = wet + 273.15;
        Tdb(n+1) = Twb(n+1);

    else if d.a==true;
            j(n + 1, 1) = DTW * fss(w(n), ima(n), Tw(n), Tdb(n), d);
            k(n + 1, 1) = DTW * gss(w(n), ima(n), Tw(n), Tdb(n), d);
            l(n + 1, 1) = DTW * hss(w(n), ima(n), Tw(n), Tdb(n), d);

            j(n + 1, 2) = DTW * fss(w(n) + (j(n + 1, 1)) / 2, ima(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
            k(n + 1, 2) = DTW * gss(w(n) + (j(n + 1, 1)) / 2, ima(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
            l(n + 1, 2) = DTW * hss(w(n) + (j(n + 1, 1)) / 2, ima(n) + (k(n + 1, 1)) / 2, Tw(n) + DTW / 2, Tdb(n), d);

            j(n + 1, 3) = DTW * fss(w(n) + (j(n + 1, 2)) / 2, ima(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
            k(n + 1, 3) = DTW * gss(w(n) + (j(n + 1, 2)) / 2, ima(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, Tdb(n), d);
            l(n + 1, 3) = DTW * hss(w(n) + (j(n + 1, 2)) / 2, ima(n) + (k(n + 1, 2)) / 2, Tw(n) + DTW / 2, Tdb(n), d);

            j(n + 1, 4) = DTW * fss(w(n) + j(n + 1, 3), ima(n) + k(n + 1, 3), Tw(n) + DTW, Tdb(n), d);
            k(n + 1, 4) = DTW * gss(w(n) + j(n + 1, 3), ima(n) + k(n + 1, 3), Tw(n) + DTW, Tdb(n), d);
            l(n + 1, 4) = DTW * hss(w(n) + j(n + 1, 3), ima(n) + k(n + 1, 3), Tw(n) + DTW, Tdb(n), d);

            w(n + 1) = w(n) + (j(n + 1, 1) + 2 * j(n + 1, 2) + 2 * j(n + 1, 3) + j(n + 1, 4)) / 6;
            ima(n + 1) = ima(n) + (k(n + 1, 1) + 2 * k(n + 1, 2) + 2 * k(n + 1, 3) + k(n + 1, 4)) / 6;
            Mep(n + 1) = Mep(n) + (l(n + 1, 1) + 2 * l(n + 1, 2) + 2 * l(n + 1, 3) + l(n + 1, 4)) / 6;

            % Calculate air temperature [K]
            Tw(n + 1) = Tw(n) + DTW;
            [~, ~, ~, ~, ~, ~, wet] = Psychrometricsnew('h',ima(n+1),'w',w(n+1));
            Twb(n+1) = wet + 273.15;
            Tdb(n+1) = Twb(n+1);
    end
    end
    %Water/Air ratio flow rate
    water(n + 1) = water(n) + d.madry * (w(n + 1) - w(n));

end


end

