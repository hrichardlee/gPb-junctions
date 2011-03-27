function out = rescale(in, p)

in = in - min(in(:));
in = in / max(in(:));

