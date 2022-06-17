
abstract type SurrogateOperatingPoint end

struct GenerationLoadScale <: SurrogateOperatingPoint
    generation_scale::Float64
    load_scale::Float64
end

function update_operating_point!(
    sys::PSY.System,
    condition::T,
) where {T <: SurrogateOperatingPoint}
    @warn "update_operating_point! not implemented for this type of surrogate operation point"
end

function update_operating_point!(sys::PSY.System, condition::GenerationLoadScale)
    generation_scale = condition.generation_scale
    load_scale = condition.generation_scale
    for g in PSY.get_components(PSY.Generator, sys)
        PSY.set_active_power!(g, PSY.get_active_power(g) * generation_scale)
    end
    for l in PSY.get_components(PSY.ElectricLoad, sys)
        PSY.set_active_power!(l, PSY.get_active_power(l) * load_scale)
        PSY.set_reactive_power!(l, PSY.get_reactive_power(l) * load_scale)
    end
end
