function i = ChooseFromLP(lp)

max_lp = max(lp);
normLogp = lp - (max_lp + log(sum(exp(lp-max_lp))));
p = exp(normLogp);
p(~isfinite(p)) = 0;
cumP = cumsum(p);
i = find(cumP>rand,1);
end

