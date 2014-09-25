function base_var = CalcBaseVar(D)
D = cast(D,'double');
off_diags = true(size(D));
for i = 1:size(D,1)
    off_diags(i,i) = false;
end
off_diags = off_diags(:);
base_var = sum((mean(D(off_diags)) - D(off_diags)).^2);

end