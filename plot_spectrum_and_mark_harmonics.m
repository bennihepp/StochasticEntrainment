function plot_spectrum_and_mark_harmonics(omega, y, om, color, options, entrainment_ratios)

min_omega = 2 * pi * options.min_frequency;
max_omega = 2 * pi * options.max_frequency;
plot_spectrum(omega, y, min_omega, max_omega);
mark_harmonics(omega, y, om, color, options, entrainment_ratios);

end
