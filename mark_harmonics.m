function mark_harmonics(omega, y, om, color, options, entrainment_ratios)

    hold on;

    S = options;

    if nargin < 4
        if isfield(S, 'entrainment_ratios')
            entrainment_ratios = S.entrainment_ratios;
        elseif isfield(S, 'entrainment_ratio')
            entrainment_ratios = S.entrainment_ratio;
        else
            entrainment_ratios = 1;
        end
    end

    for i=1:length(entrainment_ratios)
        ratio = entrainment_ratios(i);

        om_natural = 2 * pi / S.natural_period;
        om_input = om /  ratio;
        dom = S.FREQUENCY_NEIGHBOURHOOD_FACTOR * om_natural;

        for n=1:S.MAX_HARMONIC_N
            om_harmonic = om_input * n;
            [~, J] = min(abs(omega - om_harmonic));
            scatter(om_harmonic / (2*pi), abs(y(J)).^2, 100, color, 'x');
            j1 = find(omega >= om_harmonic - dom, 1, 'first');
            j2 = find(omega <= om_harmonic + dom, 1, 'last');
            if ~isempty(j1)
                scatter(omega(j1) / (2*pi), abs(y(j1)).^2, 100, color, '<');
            end
            if ~isempty(j2)
                scatter(omega(j2) / (2*pi), abs(y(j2)).^2, 100, color, '>');
            end
            if om_harmonic > 2 * pi * S.max_frequency;
                break;
            end
        end

    end

    hold off;

end
