function shape_1st_order(x)
    if x < -1
        return 0
    elseif x < 0
        return 1 + x
    elseif x < 1
        return 1 - x
    else
        return 0
    end
end

function shape_2nd_order(x)
    if x < -3/2
        return 0
    elseif x < -1/2
        return @evalpoly(x, 9/8, 3/2, 1/2)
    elseif x < 1/2
        return @evalpoly(x, 3/4, 0, -1)
    elseif x < 3/2
        return @evalpoly(x, 9/8, -3/2, 1/2)
    else
        return 0
    end
end

function shape_3rd_order(x)
    if x < -2
        return 0.
    elseif x < -1
        return @evalpoly(x, 4/3, 2, 1, 1/6)
    elseif x < 0
        return @evalpoly(x, 2/3, 0, -1, -1/2)
    elseif x < 1
        return @evalpoly(x, 2/3, 0, -1, 1/2)
    elseif x < 2
        return @evalpoly(x, 4/3, -2, 1, -1/6)
    else
        return 0.
    end
end

function shape_4th_order(x)
    if x < -5/2
        return 0.
    elseif x < -3/2
        return @evalpoly(x, 625/384, 125/48, 25/16, 5/12, 1/24)
    elseif x < -1/2
        return @evalpoly(x, 55/96, -5/24, -5/4, -5/6, -1/6)
    elseif x < 1/2
        return @evalpoly(x, 115/192, 0, -5/8, 0, 1/4)
    elseif x < 3/2
        return @evalpoly(x, 55/96, 5/24, -5/4, 5/6, -1/6)
    elseif x < 5/2
        return @evalpoly(x, 625/384, -125/48, 25/16, -5/12, 1/24)
    else
        return 0.
    end
end
