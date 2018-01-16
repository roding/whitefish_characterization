function read_input(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	output_generation_path::String = read_key(file_string, "output_generation_path", String)
	number_of_samples::Int64 = read_key(file_string, "number_of_samples", Int64)
	d::Array{Float64, 1} = read_key(file_string, "d", Array{Float64, 1})
	number_of_cells_x::Int64 = read_key(file_string, "number_of_cells_x", Int64)
	number_of_cells_y::Int64 = read_key(file_string, "number_of_cells_y", Int64)
	number_of_cells_z::Int64 = read_key(file_string, "number_of_cells_z", Int64)
	output_file_path::String = read_key(file_string, "output_file_path", String)

	return (
		output_generation_path,
		number_of_samples,
		d,
		number_of_cells_x,
		number_of_cells_y,
		number_of_cells_z,
		output_file_path)
end
