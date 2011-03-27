
%parameters that don't get sent to mex_pj_exp
pexnames = { ...
    'Thin', ...
    'Norm' ...
};
defpex = [
    1   %1  Thin
    0   %2  Norm
    ]';


pnames = {
    'Lab weight', ...
    'Positive channel weight' ...
    'Eigenvector weight', ...
    '# of slices', ...
    '# of orientations in gpb', ...
    'Radius', ...
    'Min # of angles', ...
    'Max # of angles', ...
    'Radius for thresholding', ...
    'Threshold', ...
    'Lab homogeneity weight', ...
    'Eigenvector homogeneity weight', ...
    'Hole radius', ...
    'Norm power' ...
    };
    

defp = [
    0.5    %1  Lab weight
    0.25   %2  Positive channel weight
    0.25   %3  Eigenvector weight
    16     %4  number of slices
    8      %5  number of orientations in positive channel orientation - don't change unless rerun gpb
    12     %6  radius of support
    3      %7  min number of angles
    3      %8  max number of angles
    4      %9  radius of support for calculating threshold
    0.01   %10 threshold in positive channel for pursuing further computation
    0      %11 Lab homogeneity weight
    0      %12 Eigenvector homogeneity weight
    0      %13 Inner radius of support (for hole in middle)
    1      %14 power for summing diffs/edges/homogeneities
    ]';

p.pos = @(p) [0 1 0 p(4:end)];
p.vect = @(p) [0 0 1 p(4:end)];
p.lab = @(p) [1 0 0 p(4:end)];
p.labpos = @(p) [0.5 0.5 0 p(4:end)];


p.hlab = @(p) [p(1:10) 1 0 p(13:end)];
p.hvect = @(p) [p(1:10) 0 1 p(13:end)];

p.nodiff = @(p) [0 0 0 p(4:end)];
p.nohomog = @(p) [p(1:10) 0 0];


p.wide = @(p) [p(1:5) 24 p(7:end)];
p.small = @(p) [p(1:5) 8 p(7:end)];

p.d2 = @(p) [p(1:6) 2 2 p(9:end)];

p.lowres = @(p) [p(1:3) 8 p(5:end)];

p.pow = @(x, p) [p(1:13) x];

p.hole = @(x, p) [p(1:12) x p(14:end)];



p.id = @(p) p;





p.a1 = [
    p.lab(defp)
    p.pos(defp)
    p.vect(defp)
    
    p.nodiff(p.hlab(defp))
    p.nodiff(p.hvect(defp))
    
    p.lab(p.wide(defp))
    p.pos(p.wide(defp))
    p.vect(p.wide(defp))
    
    p.pow(0.5, p.lab(defp))
    p.pow(0.5, p.pos(defp))
    p.pow(0.5, p.vect(defp))
    
    p.pow(0.25, p.lab(defp))
    p.pow(0.25, p.pos(defp))
    p.pow(0.25, p.vect(defp))
    
    p.hole(2, p.lab(defp))
    p.hole(2, p.pos(defp))
    p.hole(2, p.vect(defp))
    
    p.hole(1, p.lab(defp))
    p.hole(1, p.pos(defp))
    p.hole(1, p.vect(defp))
    ];

pex.a1 = repmat(defp, [size(p.a1, 1), 1]);
    
    
    
    
    
    