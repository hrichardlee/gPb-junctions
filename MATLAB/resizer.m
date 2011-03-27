imgsdir = 'data/2-21data/origs/';
resizeddir = 'data/resized8/';
resizefactor = 1/8;

x = dir([imgsdir, '*.JPG']);
for i=1:numel(x)
    imwrite(imresize(imread([imgsdir, x(i).name]), resizefactor), ...
        [resizeddir, x(i).name]);
end