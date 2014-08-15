ENTRAINMENT_THRESHOLD = 0.9;
MAX_HARMONIC_N = 1;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;

% VanDerPol_ArnoldTongue;

Q = zeros(size(Y, 1), size(Y, 2));
W = zeros(size(Y, 1), size(Y, 2));

for i=1:length(input_periods)
    display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
    input_period = input_periods(i);

    for j=1:length(input_amplitudes)
        display(['j=', int2str(j), ' out of ', int2str(length(input_amplitudes))]);
        input_amplitude = input_amplitudes(j);

        y = squeeze(Y(i, j, :));

        om_natural = 2 * pi / natural_period;
        om_input = 2 * pi / input_period;
        dom = FREQUENCY_NEIGHBOURHOOD_FACTOR * om_natural;

        if abs(om_natural - om_input) < 2 * dom
            Q(i, j) = inf;
        else
            % figure();
            % plot(Omega / (2*pi), abs(y).^2);
            power_total = sum(abs(y).^2);
            power_natural = compute_spectrum_power(Omega, y, om_natural, dom);
            for n=2:MAX_HARMONIC_N
                power_natural = power_natural + compute_spectrum_power(Omega, y, om_natural * n, dom);
            end
            power_input = compute_spectrum_power(Omega, y, om_input, dom);
            for n=2:MAX_HARMONIC_N
                power_input = power_input + compute_spectrum_power(Omega, y, om_input * n, dom);
            end

            Q(i, j) = power_input / power_natural;
            W(i, j) = power_input / power_total;

        end

    end

end

QQ = W';
QQ(isinf(QQ)) = max(QQ(:));

figure();
contourf(input_periods, input_amplitudes, QQ);
title('input power fraction');
xlabel('input period');
ylabel('input amplitude');
colorbar();

figure();
contourf(input_periods, input_amplitudes, log(QQ));
title('log of input power fraction');
xlabel('input period');
ylabel('input amplitude');
colorbar();

B = QQ >= ENTRAINMENT_THRESHOLD;
figure();
contourf(input_periods, input_amplitudes, B);
title(['input power fraction >= ', num2str(ENTRAINMENT_THRESHOLD)]);
xlabel('input period');
ylabel('input amplitude');
colorbar();
