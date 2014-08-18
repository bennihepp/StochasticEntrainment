PERIOD_DEVIATION_THRESHOLD = 0.01 * natural_period;
PERIODICITY_THRESHOLD = 0.05;
PERIOD_MULTIPLE_THRESHOLD = 0.05;
ENTRAINMENT_THRESHOLD = 0.9;
MAX_HARMONIC_N = 4;
FREQUENCY_NEIGHBOURHOOD_FACTOR = 0.01;

% VanDerPol_ArnoldTongue;

Q = zeros(size(Y, 1), size(Y, 2));
W = zeros(size(Y, 1), size(Y, 2));
C = zeros(size(PDmean, 1), size(PDmean, 2));


%% recompute spectras
for i=1:length(input_periods)
    display(['i=', int2str(i), ' out of ', int2str(length(input_periods))]);
    input_period = input_periods(i);

    for j=1:length(input_amplitudes)
        display(['j=', int2str(j), ' out of ', int2str(length(input_amplitudes))]);
        input_amplitude = input_amplitudes(j);

        TT = T;
        output = squeeze(X(i, j, :));

        offset_time = (tf - t0) / 5;
        offset_time = min(offset_time, 1000);
        offset = find(TT >= offset_time, 1);
        TT = TT(offset:end);
        output = output(offset:end, :);

        %% Fourier spectrum analysis
        [omega1, y1]= compute_normalized_fft_truncated(output, dt, 2*pi*min_frequency, 2*pi*max_frequency);
%         Omega = omega;
        Y(i, j, :) = y1;
    end

end


%% compute arnold tongues
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
%             for n=2:MAX_HARMONIC_N
%                 power_natural = power_natural + compute_spectrum_power(Omega, y, om_natural * n, dom);
%             end
            power_input = compute_spectrum_power(Omega, y, om_input, dom);
            power_input_harmonics = 0;
            for n=2:MAX_HARMONIC_N
                power_input_harmonics = power_input_harmonics + compute_spectrum_power(Omega, y, om_input * n, dom);
            end
            if power_input >= 0.1 * power_input_harmonics
                power_input = power_input + power_input_harmonics;
            end

            Q(i, j) = power_input / power_natural;
            W(i, j) = power_input / power_total;

        end

        mean_peak_distance = PDmean(i, j);
        std_peak_distance = PDstd(i, j);
        mean_period = mean_peak_distance * dt;
        if abs(mean_period - input_period) < PERIOD_DEVIATION_THRESHOLD
            C(i, j) = std_peak_distance / mean_peak_distance < PERIODICITY_THRESHOLD;
        end
%         factor = mean_period / input_period;
%         if mean_period < input_period
%             factor = input_period / mean_period;
%         end
%         if abs(factor - round(factor)) < PERIOD_MULTIPLE_THRESHOLD
%             C(i, j) = std_peak_distance / mean_peak_distance < PERIODICITY_THRESHOLD;
%         end

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
contourf(input_periods, input_amplitudes, B, 1);
title(['input power fraction >= ', num2str(ENTRAINMENT_THRESHOLD)]);
xlabel('input period');
ylabel('input amplitude');
colorbar();

CC = C';
figure();
contourf(input_periods, input_amplitudes, CC, 1);
title('arnold tongue');
xlabel('input period');
ylabel('input amplitude');
colorbar();
