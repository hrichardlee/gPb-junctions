%params: im, pjs, angs
%relies on subplot1 (currently in Documents/MATLAB)

function visualize(im, pjs, angs, radius)

mag = 500;

thefig = figure;
subplot1(1, 2);
subplot1(1);
imshow(im, 'InitialMagnification', mag);
subplot1(2);
imshow(pjs, [], 'InitialMagnification', mag);

xdim = size(im, 1);
ydim = size(im, 2);


key = 0;
while (key ~= 27)
	[x, y, key] = ginput(1);
	x = round(x);
	y = round(y);
	if (x > 0 && x <= xdim && y > 0 && y <= ydim)
    subplot1(1);
		imshow(im, 'InitialMagnification', mag);
		hold on;
		
		disp([num2str(x) ', ' num2str(y)]);
		
		for i = angs{y, x}
			plot([x, x - radius * sin(i)], [y, y - radius * cos(i)]);
		end
		
		hold off;
        subplot1(2);
	end
end

end
