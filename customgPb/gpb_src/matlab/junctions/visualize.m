%params: im, pjs, angs

function visualize(im, pjs, angs, radius)

close all;

imfig = figure;
imshow(im, 'InitialMagnification', 250);
jfig = figure;
imshow(pjs, 'InitialMagnification', 250);

xdim = size(im, 1);
ydim = size(im, 2);

key = 0;
while (key ~= 27)
	[x, y, key] = ginput(1);
	x = round(x);
	y = round(y);
	if (x > 0 && x <= xdim && y > 0 && y <= ydim)
		figure(imfig);
		imshow(im);
		hold on;
		
		disp([num2str(x) ', ' num2str(y)]);
		
		for i = angs{y, x}
			plot([x, x - radius * sin(i)], [y, y - radius * cos(i)]);
		end
		
		hold off;
		figure(jfig);
	end
end

end
