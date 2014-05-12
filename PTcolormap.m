function M = PTcolormap(n, range)
% Implementation of Paul Tol's divergent color scheme
% http://www.sron.nl/~pault/colourschemes.pdf

if (nargin < 2)
    steps = linspace(0,1,n);
else
    neg_prop = -1*range(1)/(range(2)-range(1));
    if (neg_prop > 0)
        steps = feval(@(x) (x < neg_prop).*(0.5/neg_prop*x) + (x >= neg_prop).*(0.5/(1-neg_prop)*(x-neg_prop) + 0.5), ...
            linspace(0,1,n));
    else
        steps = linspace(0.5,1,n);
    end
end
M = zeros(n,3);
M(:,1) = feval(@(x) 0.237 - 2.13*x + 26.92*x.^2 - 65.5*x.^3 + 63.5*x.^4 - 22.36*x.^5, steps);
M(:,2) = feval(@(x) ((0.572 + 1.524*x - 1.811*x.^2)./(1 - 0.291*x + 0.1574*x.^2)).^2, steps);
M(:,3) = feval(@(x) 1./(1.579 - 4.03*x + 12.92*x.^2 - 31.4*x.^3 + 48.6*x.^4 - 23.36*x.^5), steps);

end