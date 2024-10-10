```@meta
CurrentModule = PowerSimulationsDynamicsSurrogates
```

# PowerSimulationsDynamicsSurrogates.jl

## Background

As the name suggests, [PowerSimulationsDynamicsSurrogates.jl](https://github.com/m-bossart/PowerSimulationsDynamicsSurrogates.jl) extends [PowerSimulationsDynamics.jl](https://github.com/NREL-SIIP/PowerSimulationsDynamics.jl) for integrating surrogate models into dynamic power system simulations. On a high level, there are three major tasks assocaited with using surrogates in this application:

1) Generating training data.
2) Training the surrogate.
3) Integrating the trained surrogate into simulations.

The scope of this package is to provide generic code for tasks 1 and 3 above. Due to the diverse structure and training techniques for surrogates, code related to task 2 is left to dedicated respositories.

## Generating Training Data

### System Manipulation

When training a surrogate for part of a larger system, it can be useful to create a new `System` which consists of only part of the original system. This can be accomplished with the `create_subsystem_from_buses` function:

```@docs
create_subsystem_from_buses
```

### Generating Data

The high level function for generating data is `generate_surrogate_data`. The user specifies the operating point of the system, the types of perturbations, and the type of data to be collected:

```@docs
generate_surrogate_data
```

### Adding New Methods

There are three key areas where the inputs for `generate_surrogate_data` can be extended:

1) Adding a new perturbation:

- Define a `struct NewPerturbation <: SurrogatePerturbation`.
- Implement `add_surrogate_perturbation!` method for the new struct.

2) Adding a new operating point:

- Define a `struct NewOperatingPoint <: SurrogateOperatingPoint`.
- Implement `update_operating_point!` method for the new struct.

3) Adding a new type of dataset:

- Define a  `mutable struct NewDataset <: SurrogateDataset` .
- Implement `fill_surrogate_data!` method for the new struct.

## Integrating Trained Surrogates

`PowerSimulationsDynamicsSurrogates.jl` overloads methods from both `PowerSystems.jl` `and PowerSimulationsDynamics.jl` so that surrogate models can be added to `System` and simulated just like any other model. Refer to the documentation for details on how to add new models.

## Example

See `test/test_data_generation.jl` for an example of generating training data.

```@contents

```
