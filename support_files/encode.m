function x = encode(u)
    N = size(u, 1);
    assert(N >= 1 && mod(log2(N), 1) == 0, ...
           'The number of rows must be a power of two.');

    x = u;

    for half = 2.^(0:log2(N)-1)
        left = reshape((1:half)' + (0:2*half:N-1), [], 1);
        x(left, :) = bitxor(x(left, :), x(left + half, :));
    end
end