function plot_heatmap(x, y, values, type, levels)
    if type == 'surface'
        surf(x, y, double(values),  'LineStyle', 'none');
        xlim([min(x), max(x)]);
        ylim([min(y), max(y)]);
        view(0, 90);
    elseif type == 'contourf'
        contourf(x, y, values, levels);
    else
        error(['invalid plot type: ', type]);
    end
end
