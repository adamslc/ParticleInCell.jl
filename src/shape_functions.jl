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
