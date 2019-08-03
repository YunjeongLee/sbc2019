function penalty_term = AIC_penalty(num_data, num_par)

nd = num_data;
np = num_par;

f = @(x) 2 .* (x + 1) + 2 .* (x + 1) .* (x + 2)./(nd - x - 2);

x = 1:(nd - 3);
val = f(x);
xq = 1:20;
pf = polyfit( x, val, 2 );
fval = polyval(pf,xq);

if np <= (nd - 3)
    penalty_term = val(np);
else
    penalty_term = fval(np);
end

end