function power = compute_spectrum_power(omega, y, om, dom)

    i1 = find(omega >= om - dom, 1, 'first');
    i2 = find(omega <= om + dom, 1, 'last');

    if i1 == i2 || i1 > i2
%         warning('using only one entry for power computation');
        [~, i1] = min(abs(omega - om));
        i2 = i1;
    end

    power = sum(abs(y(i1:i2) .^ 2));

end