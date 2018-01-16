function read_output(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	number_of_samples::Int64 = read_key(file_string, "number_of_samples", Int64)
	d::Array{Float64, 1} = read_key(file_string, "d", Array{Float64, 1})
	S2::Array{Float64, 1} = read_key(file_string, "S2", Array{Float64, 1})
	t_exec::Float64 = read_key(file_string, "execution_time", Float64)

	return (
		number_of_samples,
		d,
		S2,
		t_exec)
end
