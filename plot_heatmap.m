function plot_heatmap(x, y, values, type, levels)
    if type == 'surface'
        imagesc(x, y(end:-1:1), values(end:-1:1, :));
%         surf(x, y, double(values),  'LineStyle', 'none');
%         xlim([min(x), max(x)]);
%         ylim([min(y), max(y)]);
%         view(0, 90);
    elseif type == 'contourf'
        contourf(x, y(end:-1:1), values, levels);
    elseif type == 'matrix'
        imagesc(values);
    else
        error(['invalid plot type: ', type]);
    end
end
