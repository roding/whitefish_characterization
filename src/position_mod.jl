function position_mod(x::Float64, L::Float64)
	return mod(x, L) # For characterization we need to do full modulus operation.
end
