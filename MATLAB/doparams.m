function ret = doparams(imgdir, imgpattern, gpbdir, outdir, ...
    pnames, allparams, allpexs, basenum)

if nargin < 8
    basenum = 0;
end

%so all params has (number of parameters) columns, and (number of tests)
%rows, and each row is a parameter set


addpath('~/Documents/customgPb/globalPb/lib');
addpath('~/Documents/customgPb/gpb_src/matlab/junctions');

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

interparams = {};

x = dir([imgdir, imgpattern]);
for i=1:numel(x)
    imgname = [imgdir, x(i).name];
    gpbname = [gpbdir, x(i).name, '.mat'];
    
    %do gpb unless it's already done
    if exist(gpbname, 'file')
        disp(['gpb already done for ', imgname]);
    else
        disp(['computing gpb for ', imgname]);
        [thin, gpb, mpb, spb, vect, maxo, thinmaxo] = globalPb(imgname);
        save(gpbname, 'thin', 'gpb', 'mpb', 'spb', 'vect', 'maxo', 'thinmaxo');
    end
    
    disp(['compiling pj for ', imgname]);
    for j = 1:size(allparams, 1)
        pjname = [outdir, x(i).name, '.', num2str(basenum + j), '.mat'];
        interparams{end + 1} = ...
            {imgname, gpbname, allparams(j, :), allpexs(j, :), pjname}; %params for interpj
    end
end

disp(['Starting parallel run:']);
ret = dfeval(@interpj, interparams, 'Configuration', 'local');

paramsfile = fopen([outdir, 'params.txt'], 'w');
for j = 1:size(allparams, 1)
    fprintf(paramsfile, '=====Parameters %d\n', j);
    for i = 1:size(allparams, 2)
        fprintf(paramsfile, '%s: %d\n', pnames{i}, allparams(j, i));
    end
    fprintf(paramsfile, '\n\n');
end
fclose(paramsfile);

end
