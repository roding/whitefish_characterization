function read_output(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	number_of_samples::Int64 = read_key(file_string, "number_of_samples", Int64)
	d2::Array{Float64, 1} = read_key(file_string, "d2", Array{Float64, 1})
	S2::Array{Float64, 1} = read_key(file_string, "S2", Array{Float64, 1})
	d3::Array{Float64, 1} = read_key(file_string, "d3", Array{Float64, 1})
	theta3::Array{Float64, 1} = read_key(file_string, "theta3", Array{Float64, 1})
	S3::Array{Float64, 1} = read_key(file_string, "S3", Array{Float64, 1})
	t_exec::Float64 = read_key(file_string, "execution_time", Float64)

	return (
		number_of_samples,
		d2,
		S2,
		d3,
		theta3,
		S3,
		t_exec)
end
