function showChannels(im)

nchan = size(im, 3);
cols = ceil(sqrt(nchan));
rows = ceil(nchan/cols);

for i = 1:nchan
    subplot(rows, cols, i)
    imshow(im(:,:,i))
end

end
