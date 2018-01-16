function write_input(	file_path::String,
						output_generation_path::String,
						number_of_samples::Int64,
						d::Array{Float64, 1},
						number_of_cells_x::Int64,
						number_of_cells_y::Int64,
						number_of_cells_z::Int64,
						output_file_path::String)

	file_stream::IOStream = open(file_path, "w")

	@printf(file_stream, "%s", "<input_characterization>\n")

	write_key(file_stream, "output_generation_path", output_generation_path)
	write_key(file_stream, "number_of_samples", number_of_samples)
	write_key(file_stream, "d", d)
	write_key(file_stream, "number_of_cells_x", number_of_cells_x)
	write_key(file_stream, "number_of_cells_y", number_of_cells_y)
	write_key(file_stream, "number_of_cells_z", number_of_cells_z)
	write_key(file_stream, "boundary_condition", boundary_condition)
	write_key(file_stream, "output_file_path", output_file_path)

	@printf(file_stream, "%s", "</input_characterization>")

	close(file_stream)

	nothing
end
