"""
    Qvariable <: IndependentVariable

A variable tied to the force density of an FDM element

```julia
QVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
QVariable(element::Asap.AbstractElement, value::Float64, lowerbound::Float64, upperbound::Float64)
QVariable(element::Asap.AbstractElement, lowerbound::Float64, upperbound::Float64)
```
"""
mutable struct QVariable <: IndependentVariable
    i::Int64
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Int64

    function QVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64)
        new(elementindex, value, lowerbound, upperbound, 0)
    end

    function QVariable(element::FDMelement, value::Float64, lowerbound::Float64, upperbound::Float64)
        new(element.elementID, value, lowerbound, upperbound, 0)
    end

    function QVariable(element::FDMelement, lowerbound::Float64, upperbound::Float64)
        new(element.elementID, value, lowerbound, upperbound, 0)
    end
end