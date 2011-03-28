
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
    'L weight', ...                         1
    'ab weight', ...                        2
    'Positive channel weight' ...           3
    'Eigenvector weight', ...               4
    'Texture weight', ...                   5
    ...
    'L homogeneity weight', ...             6
    'ab homogeneity weight', ...            7
    'Eigenvector homogeneity weight', ...   8
    'Texture homgeneity weight', ...        9
    ...
    '# of slices', ...                      10
    '# of orientations in gpb', ...         11
    ...
    'Radius1', ...                          12
    'Radius2', ...                          13
    'Hole radius1', ...                     14
    ...
    'Min # of angles', ...                  15
    'Max # of angles', ...                  16
    ...
    'Radius for thresholding', ...          17
    'Threshold', ...                        18
    ...
    'Norm power' ...                        19
    };
    

defp = [
    1      %1  L weight
    0      %2  ab weight
    0      %3  Positive channel weight
    0      %4  Eigenvector weight
    0      %5  Texture weight
    
    0      %6  L homog weight
    0      %7  ab homog weight
    0      %8  eigenvector homog weight
    0      %9  texture homog weight
    
    16     %10 number of slices
    8      %11 number of orientations in positive channel orientation - don't change unless rerun gpb
    
    12     %12 radius of support (1)
    0      %13 radius of support (2)
    0      %14 hole radius
    
    3      %15 min number of angles
    3      %16 max number of angles
    
    4      %17 radius of support for calculating threshold
    0.01   %18 threshold in positive channel for pursuing further computation
 
    1      %19 power for summing diffs/edges/homogeneities
    ]';

p.pos = @(p)    [0 0 1 0 0 0 0 0 0 p(10:end)];
p.vect = @(p)   [0 0 0 1 0 0 0 0 0 p(10:end)];
p.lab = @(p)    [0.5 0.5 0 0 0 0 0 0 0 p(10:end)];
p.l = @(p)      [1 0 0 0 0 0 0 0 0 p(10:end)];
p.ab = @(p)     [0 1 0 0 0 0 0 0 0 p(10:end)];
p.text = @(p)   [0 0 0 0 1 0 0 0 0 p(10:end)];

p.hlab = @(p)   [0 0 0 0 0 0.5 0.5 0 0 p(10:end)];
p.hl = @(p)     [0 0 0 0 0 1 0 0 0 p(10:end)];
p.hab = @(p)    [0 0 0 0 0 0 1 0 0 p(10:end)];
p.hvect = @(p)  [0 0 0 0 0 0 0 1 0 p(10:end)];

p.wide = @(p) [p(1:11) 24 p(13:end)];
p.small = @(p) [p(1:11) 8 p(13:end)];
p.wide2 = @(p) [p(1:12) 24 p(14:end)];

p.hole = @(x, p) [p(1:13) x p(15:end)];

p.d2 = @(p) [p(1:14) 2 2 p(17:end)];

p.pow = @(x, p) [p(1:18) x];





p.a1 = [
    p.lab(defp)
    p.pos(defp)
    p.vect(defp)
    
    p.hlab(defp)
    p.hvect(defp)
    
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
    
p.a2 = [p.a1(1:5, :); p.a1(9:end, :)];
for i=1:size(p.a2, 1)
    p.a2(i, :) = p.wide2(p.a2(i, :));
end
pex.a2 = repmat(defp, [size(p.a2, 1), 1]);

    
    
    