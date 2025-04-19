function out = mmod(a, b)

    %MMOD   Minmod slope limiter (vectorized)
    out = sign(a) .* max(0, min(abs(a), sign(a).*b));

end