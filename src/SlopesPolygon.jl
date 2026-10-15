
"""
Abstract type that describes a Newton Polygon/Hodge Polygon
like object.

There are a couple of ways to store such data, the idea
is that different implementations can be used to implement
those different strategies
"""
abstract type AbstractPolygon end

"""
    SlopesPolygon(slopelengths)
    SlopesPolygon(xs, yvals)
    SlopesPolygon(vertices)

Represent a piecewise-linear polygon starting at `(0, 0)`. The integer-vector
constructor uses slopes `0, 1, 2, ...`; each entry of `slopelengths` is the
horizontal length for the corresponding slope. For example,
`SlopesPolygon([1, 19, 1])` is the Hodge polygon of the primitive second
cohomology of a K3 surface.

The two-vector constructor takes points `(xs[i], yvals[i])` with `xs` sorted
in ascending order and constructs their lower convex hull. The vertices
constructor takes the polygon's vertices, including the initial `(0, 0)`.

The object stores its slopes, their lengths, the values at integer
coordinates, and the slope immediately before each coordinate. Indexing with
`sp[x]` returns the polygon's value at the integer coordinate `x`.

Returns a `SlopesPolygon` value. Its representation is computed eagerly.

### Examples

```jldoctest
julia> polygon = DeRham.SlopesPolygon([1, 19, 1]);

julia> polygon[2] == 1
true
```
"""
struct SlopesPolygon
    slopes::Array{Rational{Int}}
    slopelengths::Array{Int}
    values::Array{Rational{Int}}
    slopesbefore::Array{Rational{Int}}
end

# MARK - constructors

"""
Gives the array of values of the polygon,
starting with zero. Thus,
(i-1,values(sp)[i])
is a point on the polygon.

Right now, we call this in the constructors.
"""
function values(slopes, slopelengths)
    l = sum(slopelengths) + 1
    vals = zeros(Rational{Int}, l)
    slopesbefore = zeros(Rational{Int}, l)

    vals[1] = 0 # starting point is (0,0)
    slopesbefore[1] = 0

    nSlopes = length(slopelengths)
    j = 1
    for i = 1:nSlopes
        slope = slopes[i]
        slopeend = j-1 + slopelengths[i]
        while j <= slopeend
            yval = vals[j] + slope
            vals[j+1] = yval
            slopesbefore[j+1] = slope
            j = j + 1
        end
    end

    (vals, slopesbefore)
end

"""
    SlopesPolygon(vertices::Vector{Tuple{Int64,Int64}})

Construct a polygon from its vertices, beginning with `(0, 0)`. Each
successive pair determines a slope and horizontal length.

- `vertices`: vertices in increasing `x`-coordinate order.

Returns a `SlopesPolygon`.
"""
function SlopesPolygon(vertices::Vector{Tuple{Int64,Int64}})
    #TODO implement

    n = length(vertices)-1
    slopelengths = zeros(Int, n)
    slopes = zeros(Rational{Int}, n)
    for i = 2:(n+1)
        slopevec = vertices[i] .- vertices[i-1]
        slopes[i-1] = slopevec[2] // slopevec[1]
        slopelengths[i-1] = slopevec[1]
    end

    SlopesPolygon(slopes, slopelengths, values(slopes, slopelengths)...)
end

"""
    SlopesPolygon(slopelengths::Vector{<:Integer})

Construct a polygon with slopes `0, 1, 2, ...`, where `slopelengths[i]` is
the horizontal length of the `i`-th slope. For example,
`SlopesPolygon([1, 19, 1])` is the Hodge polygon of the primitive second
cohomology of a K3 surface.

- `slopelengths`: nonnegative integer lengths, one for each slope.

Returns a `SlopesPolygon`.
"""
function SlopesPolygon(slopelengths::Array{T}) where {T<:Integer}
    n = length(slopelengths) - 1
    slopes = Rational.(collect(0:n))
    SlopesPolygon(slopes, slopelengths, values(slopes, slopelengths)...)
end

"""
    SlopesPolygon(xs::Vector{T}, yvals::Vector{T}) where {T<:Real}

Construct the lower convex hull of the points `(xs[i], yvals[i])`. The input
assumes the first point is `(0, 0)` and the x-coordinates are sorted.

- `xs`: x-coordinates, sorted in ascending order.
- `yvals`: corresponding y-coordinates.

Returns a `SlopesPolygon` whose vertices lie on the lower convex hull.
"""
function SlopesPolygon(xs::Array{T}, yvals::Array{T}) where {T<:Real}
    vertices = [(0, 0)]

    points = [(xs[i], yvals[i]) for i in eachindex(yvals)]

    x = 0
    while x < xs[end] # i.e. end of the hodge polygon

        startpoint = points[x+1]
        x2 = x + 1
        smallest_slope = typemax(Int128)
        finalx = x

        while x2 ≤ length(yvals) - 1

            thispoint = points[x2+1]
            thisslope = (thispoint[2] - startpoint[2]) // (thispoint[1] - startpoint[1])
            if thisslope ≤ smallest_slope # we do \leq even though the slope won't change so we can update finalpoint
                smallest_slope = thisslope
                finalx = x2
            end
            x2 = x2 + 1
        end # postcondition: x2 is the next vertex

        # add x2 to the vertices
        push!(vertices, points[finalx+1])

        x = finalx
    end

    SlopesPolygon(vertices)
end

function SlopesPolygon(coefficients, valuation)
    n = length(coefficients)-1
    SlopesPolygon(collect(0:n), reverse!(valuation.(coefficients)))
end

# MARK - creating polygons from other polygons

tatetwist(sp::SlopesPolygon, n) = SlopesPolygon(sp.slopes .- n, sp.slopelengths, nothing)

# MARK - methods

slopelengths(sp::SlopesPolygon) = sp.slopelengths

function Base.getindex(sp::SlopesPolygon, i::Integer)
    sp.values[i+1]
end

function Base.:(==)(x::SlopesPolygon, y::SlopesPolygon)
    (x.slopes == y.slopes) && (x.slopelengths == y.slopelengths)
end

function Base.hash(x::SlopesPolygon, h::UInt64)
    hash((x.slopes, x.slopelengths))
end

endcoord(sp::SlopesPolygon) = length(values(sp)) - 1

function vertices(sp::SlopesPolygon)
    error("not implemented")
end
