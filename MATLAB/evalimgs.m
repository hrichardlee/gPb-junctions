function evalimgs(matdir, imgdir, radius, threshold, vectorout)

if nargin == 4
    vectorout = true;
end

if vectorout
    outtype = 'pdf';
    outcommand = '-dpdf';
else
    outtype = 'bmp';
    outcommand = '-dbmp';
end

x = dir([matdir, '*.mat']);
for i_img = 1:numel(x)
    %%%make names
    shortname = x(i_img).name;
    matname = [matdir, shortname];
    dotindex = strfind(shortname, '.');
    origname = [imgdir, shortname(1:(dotindex(2) - 1))];
    
    load(matname); %this gives us pjs, angs, ests, params, pex
    %string that summarizes parameters
    pstr = [params2text(params, pex)  ...
        '-Nmr' num2str(radius) 'th' num2str(threshold)];
    
    outname = [matdir, shortname(1:end - 3), pstr, '.', outtype];
    
    %get image
    im = im2double(imread(origname));
    
    %find nonmaxes
    [r, c] = nonmaxsuppts(pjs, radius, threshold);
    
    %find range of nonmax pjs
    nonmaxinds = sub2ind(size(pjs), r, c);
    nonmaxpjs = pjs(nonmaxinds);
    minpj = min(nonmaxpjs);
    maxpj = max(nonmaxpjs);
    
    %add min max range of pjs to the title, but not the filename
    pstr = [pstr, '-', num2str(minpj, 2), 'to', num2str(maxpj, 2)];
    
    %draw nonmaxes
    figure(1), imshow(im), hold on
    title(pstr);
    for i_pt = 1:numel(c)
        ri = r(i_pt);
        ci = c(i_pt);
    	for i_ang = angs{ri, ci}
            currpj = pjs(ri, ci);
            currrad = 2 + 20 * (currpj - minpj) / (maxpj - minpj);
    		plot([ci, ci - currrad * sin(i_ang)], [ri, ri - currrad * cos(i_ang)]);
    	end
    end
    hold off
    print(1, outname, outcommand);
end


end



%%%%just finds points of nonmaxes
function [r,c] = nonmaxsuppts(cim, radius, thresh)
    
    % Extract local maxima by performing a grey scale morphological
    % dilation and then finding points in the corner strength image that
    % match the dilated image and are also greater than the threshold.
    
    sze = 2*radius+1;                   % Size of dilation mask.
    mx = ordfilt2(cim,sze^2,ones(sze)); % Grey-scale dilate.

    % Make mask to exclude points within radius of the image boundary. 
    bordermask = zeros(size(cim));
    bordermask(radius+1:end-radius, radius+1:end-radius) = 1;
    
    % Find maxima, threshold, and apply bordermask
    cimmx = (cim==mx) & (cim>thresh) & bordermask;
    
    [r,c] = find(cimmx);                % Find row,col coords.
    
end

    
