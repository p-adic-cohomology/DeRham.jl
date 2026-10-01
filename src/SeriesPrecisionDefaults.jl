# SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

# Table values from controlledreduction (https://github.com/edgarcosta/controlledreduction, tools/default_args.cc), Copyright (C) 2013-2017 Edgar Costa.

# Copyright (c) 2024 the DeRham.jl contributors.

"""
    series_precision_default(p, n, d)

Tabulated series-precision bounds for small `(n, d)`, taken from
`default_args()` in controlledreduction's `tools/default_args.cc`.
Returns the `N` vector for the given prime `p`, or `nothing` when `(n, d, p)`
has no table entry.
"""
function series_precision_default(p, n, d)
    N = nothing

    if (n == 2 && d == 3)
        if (p < 5)
            N = [2, 3]
        elseif (5 <= p && p < 17)
            N = [2, 2]
        elseif (17 <= p)
            N = [1, 1]
        end
    elseif (n == 2 && d == 4)
        if (p < 5)
            N = [4, 4]
        elseif (5 <= p && p < 17)
            N = [3, 3]
        elseif (17 <= p)
            N = [2, 2]
        end
    elseif (n == 2 && d == 5)
        if (p < 5)
            N = [6, 6]
        elseif (5 <= p && p < 7)
            N = [4, 5]
        elseif (7 <= p)
            N = [4, 4]
        end
    elseif (n == 2 && d == 6)
        if (p < 5)
            N = [8, 9]
        elseif (5 <= p && p < 7)
            N = [7, 7]
        elseif (7 <= p && p < 11)
            N = [6, 7]
        elseif (11 <= p)
            N = [6, 6]
        end
    elseif (n == 2 && d == 7)
        if (p < 5)
            N = [11, 11]
        elseif (5 <= p && p < 7)
            N = [10, 10]
        elseif (7 <= p && p < 11)
            N = [10, 10]
        elseif (11 <= p && p < 17)
            N = [9, 9]
        elseif (17 <= p)
            N = [8, 8]
        end
    elseif (n == 2 && d == 8)
        if (p < 5)
            N = [14, 14]
        elseif (5 <= p && p < 7)
            N = [13, 13]
        elseif (7 <= p && p < 13)
            N = [13, 13]
        elseif (13 <= p && p < 17)
            N = [12, 13]
        elseif (17 <= p)
            N = [11, 11]
        end
    elseif (n == 2 && d == 9)
        if (p < 5)
            N = [18, 18]
        elseif (5 <= p && p < 7)
            N = [16, 16]
        elseif (7 <= p && p < 11)
            N = [16, 16]
        elseif (11 <= p && p < 17)
            N = [16, 16]
        elseif (17 <= p)
            N = [15, 15]
        end
    elseif (n == 3 && d == 4)
        if (p < 5)
            N = [7, 7, 8]
        elseif (5 <= p && p < 7)
            N = [4, 5, 5]
        elseif (7 <= p && p < 23)
            N = [4, 4, 3]
        elseif (23 <= p && p < 43)
            N = [3, 3, 3]
        elseif (43 <= p)
            N = [3, 3, 2]
        end
    elseif (n == 3 && d == 5)
        if (p < 5)
            N = [11, 11, 10]
        elseif (5 <= p && p < 7)
            N = [8, 8, 9]
        elseif (7 <= p && p < 11)
            N = [8, 8, 7]
        elseif (11 <= p && p < 23)
            N = [7, 7, 6]
        elseif (23 <= p && p < 29)
            N = [6, 6, 6]
        elseif (29 <= p)
            N = [6, 6, 5]
        end
    elseif (n == 3 && d == 6)
        if (p < 5)
            N = [17, 18, 17]
        elseif (5 <= p && p < 7)
            N = [15, 15, 14]
        elseif (7 <= p && p < 11)
            N = [15, 15, 14]
        elseif (11 <= p && p < 17)
            N = [14, 14, 13]
        elseif (17 <= p && p < 23)
            N = [13, 13, 12]
        elseif (23 <= p)
            N = [12, 12, 11]
        end
    end

    N
end
