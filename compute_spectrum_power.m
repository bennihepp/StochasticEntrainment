function power = compute_spectrum_power(omega, y, om, dom)

    i1 = find(omega >= om - dom, 1, 'first');
    i2 = find(omega <= om + dom, 1, 'last');

    if isempty(i1) || isempty(i2)

        power = 0;

    else

        if i1 == i2 || i1 > i2
    %         warning('using only one entry for power computation');
            [~, I] = min(abs(omega - om));
            if isempty(I)
                power = 0;
                return;
            else
                i1 = I(1);
                i2 = I(1);
            end
        end

        power = sum(abs(y(i1:i2) .^ 2));

    end

end
