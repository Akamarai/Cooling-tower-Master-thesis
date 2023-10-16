function ima = checkima(w, Ta, d)

    [cpa, cpv, ~, ] = heatcap(Ta);
    ima = cpa * (Ta - 273.15) + w * (d.ifgwo + cpv * (Ta - 273.15));
end
