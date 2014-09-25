function r = CorrNoCentering(v1, v2)

v1 = v1(:);
v2 = v2(:);

r = v1'*v2 / (sqrt(v1'*v1) * sqrt(v2'*v2));

end

