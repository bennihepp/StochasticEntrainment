function plot_heatmap(x, y, values, type, levels)
    if strcmpi(type, 'surface')
        surf(x, y, double(values),  'LineStyle', 'none');
        xlim([min(x), max(x)]);
        ylim([min(y), max(y)]);
        view(0, 90);
    elseif strcmpi(type, 'contourf')
        contourf(x, y(end:-1:1), values, levels);
    elseif strcmpi(type, 'matrix')
        imagesc(x, y(end:-1:1), values(end:-1:1, :));
        set(gca,'YDir','normal');
    else
        error(['invalid plot type: ', type]);
    end
end
