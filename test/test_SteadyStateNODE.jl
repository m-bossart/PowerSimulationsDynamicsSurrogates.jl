@testset "SteadyStateNODE" begin
    x_scale = [1.1,1.0]
    x_bias =  [0.0, 1.0]
    initializer = Chain((x)-> x.* x_scale .+ x_bias, Dense(2,2, tanh; bias =true))
    node= Chain(Parallel(+,Dense(2,2; bias =true),Dense(2,2; bias =true)), Dense(2,2; bias =true))
    observer = Dense(2,2; bias =true)

     ssnode = SteadyStateNODE(name = "test", 
        initializer_structure=[(2, 2, true, "tanh")],
        initializer_parameters= Flux.destructure(initializer)[1],
        node_structure_exogenous=[(2, 2, true, "tanh")],
        node_structure_states=[(2, 2, true, "tanh")],
        node_structure_common=[(2, 2, true, "tanh")],
        node_parameters= Flux.destructure(node)[1],
        observer_structure=[(2, 2, true, "tanh")],
        observer_parameters= Flux.destructure(observer)[1],
        x_scale = x_scale,
        x_bias = x_bias,
        exogenous_scale =  [1.0,1.0],
        exogenous_bias =   [0.0,0.0],
    )
 
    x = rand(2)
    @test isapprox(initializer(x),initializer_psids(ssnode, x); atol = 1e-14)

end

#TODO - what automatic generation should we adopt (parse from JSON, etc.)
#TODO - where would you check that the parameters provided are the correct length for the structure provided? 
#TODO - write a test that builds a system in PSY with this component. 
#TODO - do we build the matricies every single function call to this component? 
