

sw1 = im2double(imread('swatch2.bmp'));

big = im2double(imread('30-150.bmp'));

edge = im2double(imread('edge.bmp'));


send = big;

[pjs, angs] = mex_pj_exp(send(:, :, 1), send(:, :, 2), send(:, :, 3), thin);

%js = (cellfun(@(x) numel(x), angs) == 3) .* pjs;

%degs = cellfun(@(x) rad2deg(x), angs, 'UniformOutput', false);
