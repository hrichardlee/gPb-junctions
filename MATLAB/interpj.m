
function ret = interpj(imgname, gpbname, params, pex, pjname)
    addpath('~/Documents/customgPb/globalPb/lib');
    addpath('~/Documents/customgPb/gpb_src/matlab/junctions');

    load(gpbname); %gives us gpb, thin, maxo, vect, etc
    img = im2double(imread(imgname));
    
    posimg = (pex(1) * thin) + (pex(2) * gpb);
    
    [pjs, angs, ests] = mex_pj_exp(img, posimg, maxo, vect, params);
    save(pjname, 'pjs', 'angs', 'ests', 'params', 'pex');
    
    ret = 5;
end