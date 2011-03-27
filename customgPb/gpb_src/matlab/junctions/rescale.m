function out = rescale(in, gamma)

out = in;
nonzeros = out(out ~= 0);
out = (out ~= 0) .* (out - min(nonzeros));
out = out / max(out(:));

out = out .^ gamma;

