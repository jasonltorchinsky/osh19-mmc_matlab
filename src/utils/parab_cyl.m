function res = parab_cyl(y, order)

% Returns normalized parabolic cylinder function

res = (1.0 / ((2 * pi)^(-1/4) * sqrt(2^order * factorial(order)))) ...
    * exp(-y.^2/4) .* hermiteH(order, y./sqrt(2));

end

