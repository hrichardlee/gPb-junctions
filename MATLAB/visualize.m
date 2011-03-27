%params: im, pjs, angs
%relies on subplot1 (currently in Documents/MATLAB)

function visualize(im, pjs, angs, radius, gamma)

mag = 200;

thefig = figure;
subplot1(1, 2);
subplot1(1);
imshow(im, 'InitialMagnification', mag);
subplot1(2);
imshow(rescale(pjs, gamma), 'InitialMagnification', mag);

ydim = size(im, 1);
xdim = size(im, 2);


key = 0;
while (key ~= 27)
	[x, y, key] = ginput(1);
    
    if (key == 1)
        x = round(x);
        y = round(y);
        if (x > 0 && x <= xdim && y > 0 && y <= ydim)
            subplot1(1);
            imshow(im, 'InitialMagnification', mag);
            hold on;

            disp([num2str(x) ', ' num2str(y) ': ' num2str(pjs(y, x))]);

            for i = angs{y, x}
                p = plot([x, x - radius * sin(i)], [y, y - radius * cos(i)]);
                set(p, 'Color', 'red');
            end
            
            disp(angs{y, x});

            hold off;
            subplot1(2);
        end
    elseif (key == 122)
        gamma = gamma / 1.5;
        subplot1(2);
        imshow(rescale(pjs, gamma), 'InitialMagnification', mag);
        disp(['gamma: ', num2str(gamma)]);
    elseif (key == 120)
        gamma = gamma * 1.5;
        subplot1(2);
        imshow(rescale(pjs, gamma), 'InitialMagnification', mag);
        disp(['gamma: ', num2str(gamma)]);
    else
        disp(num2str(key));
    end
end

end
