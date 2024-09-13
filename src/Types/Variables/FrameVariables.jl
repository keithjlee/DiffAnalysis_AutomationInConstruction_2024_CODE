abstract type SectionVariableProperty end
struct SectionA <: SectionVariableProperty end
struct SectionIx <: SectionVariableProperty end
struct SectionIy <: SectionVariableProperty end
struct SectionJ <: SectionVariableProperty end

"""
    SectionVariable <: IndependentVariable

A variable tied to a geometric section property of an element. 

```
SectionVariable(elementindex::Int64, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
SectionVariable(element::Element, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
SectionVariable(element::Element, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)
```
"""
mutable struct SectionVariable{T<:SectionVariableProperty} <: IndependentVariable
    i::Int64
    val::Float64
    lb::Float64
    ub::Float64
    iglobal::Float64
end

const property_to_property_type = Dict(
    :A => SectionA,
    :Ix => SectionIx,
    :Iy => SectionIy,
    :J => SectionJ
)

function SectionVariable(element::Asap.Element, value::Float64, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)

    T = property_to_property_type[property]
    
    return SectionVariable{T}(element.elementID, value, lowerbound, upperbound, 0)

end

function SectionVariable(element::Asap.Element, lowerbound::Float64, upperbound::Float64, property::Symbol = :A)

    T = property_to_property_type[property]
    
    value = getproperty(element.section, property)
    @assert lowerbound ≤ value ≤ upperbound

    return SectionVariable{T}(element.elementID, value, lowerbound, upperbound, 0)

end