function perimeter = compute_perimeter_2d(mask, dx, dy)
%COMPUTE_PERIMETER_2D  Estimate perimeter of a 2D binary mask.

    mask = logical(mask);
    perimeter = 0;
    [rows, cols] = find(mask);
    for idx = 1:length(rows)
        r = rows(idx); c = cols(idx);
        if r == 1 || ~mask(r-1, c), perimeter = perimeter + dx; end
        if r == size(mask,1) || ~mask(r+1, c), perimeter = perimeter + dx; end
        if c == 1 || ~mask(r, c-1), perimeter = perimeter + dy; end
        if c == size(mask,2) || ~mask(r, c+1), perimeter = perimeter + dy; end
    end
end
