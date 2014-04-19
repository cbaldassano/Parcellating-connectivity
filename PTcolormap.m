function M = PTcolormap(n)
% Implementation of Paul Tol's divergent color scheme
% http://www.sron.nl/~pault/colourschemes.pdf

M = zeros(n,3);
M(:,1) = feval(@(x) 0.237 - 2.13*x + 26.92*x.^2 - 65.5*x.^3 + 63.5*x.^4 - 22.36*x.^5,...
               linspace(0,1,n));
M(:,2) = feval(@(x) ((0.572 + 1.524*x - 1.811*x.^2)./(1 - 0.291*x + 0.1574*x.^2)).^2,...
               linspace(0,1,n));
M(:,3) = feval(@(x) 1./(1.579 - 4.03*x + 12.92*x.^2 - 31.4*x.^3 + 48.6*x.^4 - 23.36*x.^5),...
               linspace(0,1,n));
           
end