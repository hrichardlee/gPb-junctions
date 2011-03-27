function s = params2text(params, pex)

s = '';
s = f(s, params(1), 'Lab');
s = f(s, params(2), 'Pos');
if params(2) > 0
    s = f(s, pex(1), 'Thin');
    s = f(s, pex(2), 'Gpb');
end
s = f(s, params(3), 'Eig');
s = f(s, params(11), 'Labh');
s = f(s, params(12), 'Eigh');


s = f(s, params(6), 'Rad', 0, -1);
s = f(s, params(13), 'Holerad', 0, -1);


s = f(s, params(14), 'Pow', 1);


end


function s = f(s, param, pname, default, complete)

if nargin < 4
    default = 0;
end
if nargin < 5
    complete = 1;
end


if param ~= default
    if param == complete
        s = [s pname];
    else
        s = [s pname num2str(param)];
    end
end
   
end