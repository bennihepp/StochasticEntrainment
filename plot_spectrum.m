function plot_spectrum(omega, y, min_omega, max_omega)

j1 = find(omega >= min_omega, 1, 'first');
j2 = find(omega <= max_omega, 1, 'last');
plot(omega(j1:j2) ./ (2 * pi), y(j1:j2));
xlabel('frequency f');
ylabel('power |y|^2');

end
