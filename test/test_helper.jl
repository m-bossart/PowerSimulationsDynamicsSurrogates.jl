m = DataDrivenParams()
surr = build_data_driven_model(m, SteadyStateNODE, "test") 
surr2 = build_data_driven_model(m, SteadyStateNODEObs, "test")




standard_inverter =  DynamicInverter(
    name = "test-name",
    ω_ref = 1.0, 
    converter = converter_high_power(), 
    outer_control = outer_control(), 
    inner_control = inner_control(), 
    dc_source = dc_source_lv(), 
    freq_estimator = pll(), 
    filter = filt(), 
)

dd_pll_params = DataDrivenParams(
    input_normalization = true, 
    target_normalization = true, 
    input_ref_frame = true, 
    target_ref_frame = true, 
    ) 

dd_inverter =  DynamicInverter(
    name = "test-name",
    ω_ref = 1.0, 
    converter = converter_high_power(), 
    outer_control = outer_control(), 
    inner_control = inner_control(), 
    dc_source = dc_source_lv(), 
    freq_estimator = build_data_driven_model(DataDrivenParams(), DataDrivenFrequencyEstimator, "test-pll"),
    filter = filt(), 
)

#We might need to specify if we wan

build_data_driven_model