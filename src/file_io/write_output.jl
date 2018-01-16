function write_output(	file_path::String,
						number_of_samples::Int64,
						d::Array{Float64, 1},
						S2::Array{Float64, 1},
						t_exec::Float64)

	file_stream::IOStream = open(file_path, "w")

	@printf(file_stream, "%s", "<output_characterization>\n")

	write_key(file_stream, "number_of_samples", number_of_samples)
	write_key(file_stream, "d", d)
	write_key(file_stream, "S2", S2)
	write_key(file_stream, "execution_time", t_exec)

	@printf(file_stream, "%s", "</output_characterization>")

	close(file_stream)

	nothing
end
