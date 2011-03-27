function ret = doimgs(imgdir, imgpattern, gpbdir, pjdir, params)

addpath('~/Documents/customgPb/globalPb/lib');
addpath('~/Documents/customgPb/gpb_src/matlab/junctions');

logfile = fopen('log.log', 'a');

%imgdir = 'models/juncsorig/';
%gpbdir = 'models/juncsgpb/';
%pjdir = 'models/juncspj/';
%imgpattern = '*.bmp';

x = dir([imgdir, imgpattern]);
for i=1:numel(x)
    imgname = [imgdir, x(i).name];
    gpbname = [gpbdir, x(i).name, '.mat'];
    pjname = [pjdir, x(i).name, '.mat'];
    pjimgname = [pjdir, x(i).name, '.bmp'];
    
    %do gpb unless it's already done
    if exist(gpbname, 'file')
        disp(['gpb already done for ', imgname]);
        tgpb = 0;
        load(gpbname);
    else
        disp(['computing gpb for ', imgname]);
        gpbstartt = tic;
        [thin, gpb, mpb, spb, vect, maxo, thinmaxo] = globalPb(imgname);
        tgpb = toc(gpbstartt);
        save(gpbname, 'thin', 'gpb', 'mpb', 'spb', 'vect', 'maxo', 'thinmaxo');
    end
    
    %do pj unless it's already done
    %(so we can continue interrupted computations with granularity
    %of one image)
    if exist(pjname, 'file')
        disp(['pj already done for ', imgname]);
        tpj = 0;
    else
        disp(['computing pj for ', imgname]);
        img = im2double(imread(imgname));

        pjstartt = tic;
        [pjs, angs, ests] = mex_pj_exp(img, thin, thinmaxo, vect, params);
        tpj = toc(pjstartt);
        disp(['time for pj: ', num2str(tpj)]);
        save(pjname, 'pjs', 'angs', 'ests', 'params');
        
        normpjs = pjs / max(max(pjs));
        imwrite(normpjs, pjimgname);
    end
    
    %write log
    fprintf(logfile, '%s: %d %d\n', x(i).name, tgpb, tpj);
end

fclose(logfile);

ret = 5;