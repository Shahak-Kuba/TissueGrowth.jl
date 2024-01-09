using Base

function Base.insert!(A::ElasticMatrix, col, xy)
    i1 = 2col - 1
    i2 = 2col
    data = A.data
    insert!(data, i1, xy[1])
    insert!(data, i2, xy[2])
    return A
end

function Base.deleteat!(A::ElasticMatrix, col)
    data = A.data
    deleteat!(data,2col-1)
    deleteat!(data,2col-1)
    return A
end