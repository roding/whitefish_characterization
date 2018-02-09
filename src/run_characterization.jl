workspace()

include("file_io/read_key.jl")
include("file_io/write_key.jl")

include("file_io/read_input.jl")
include("file_io/read_output_generation.jl")
include("file_io/write_output.jl")

include("characteristic_matrix_ellipse.jl")
include("characteristic_matrix_ellipsoid.jl")
include("inverse_characteristic_matrix_ellipsoid.jl")
include("rotation_matrix.jl")
include("inverse_rotation_matrix.jl")

include("generate_cell_lists.jl")
include("overlap_cuboid_binary.jl")
include("axis_aligned_bounding_box.jl")
include("intersect_box_box.jl")

foo = @__FILE__
@eval @everywhere f = $foo
#@everywhere println(f)
@everywhere (program_file_dir, program_file_name) = splitdir(f)
@everywhere include(joinpath(program_file_dir, "characterize_S2.jl"))
@everywhere include(joinpath(program_file_dir, "characterize_S3.jl"))
@everywhere include(joinpath(program_file_dir, "characterize_zeta2.jl"))
@everywhere include(joinpath(program_file_dir, "signed_distance_mod.jl"))
@everywhere include(joinpath(program_file_dir, "position_mod.jl"))
@everywhere include(joinpath(program_file_dir, "rotation_matrix.jl"))
@everywhere include(joinpath(program_file_dir, "generate_random_unit_quaternion.jl"))
@everywhere include(joinpath(program_file_dir, "quaternion_mult.jl"))
@everywhere include(joinpath(program_file_dir, "rotate.jl"))

function run_characterization()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Start time.
	t_start_ns::Int64 = convert(Int64, time_ns())

	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)

	# Assert that input is file and store path.
	input_file_path::String = ""
	if isfile(ARGS[1])
		input_file_path = ARGS[1]
	else
		println("No input file specified or specified input file does not exist. Aborting.")
		return nothing
	end

	# Read characerization input from file.
	println(join(("Reading characerization input from file ", input_file_path, "...")))
	(	output_generation_path::String,
		number_of_samples::Int64,
		d2::Array{Float64, 1},
		d3::Array{Float64, 1},
		theta3::Array{Float64, 1},
		number_of_cells_x::Int64,
		number_of_cells_y::Int64,
		number_of_cells_z::Int64,
		output_file_path::String) = read_input(input_file_path)

	# Read generation output from file.
	println(join(("Reading generation output from file ", output_generation_path, "...")))
	(	particle_type::String,
		R::Array{Float64, 2},
		Lx::Float64,
		Ly::Float64,
		Lz::Float64,
		phi::Float64,
		X::Array{Float64, 1},
		Y::Array{Float64, 1},
		Z::Array{Float64, 1},
		Q0::Array{Float64, 1},
		Q1::Array{Float64, 1},
		Q2::Array{Float64, 1},
		Q3::Array{Float64, 1},
		execution_time_generation::Float64) = read_output_generation(output_generation_path)
	number_of_particles::Int64 = length(X)

	# Characteristic/rotation matrix entries.
	A11::Array{Float64, 1} = zeros(number_of_particles)
	A12::Array{Float64, 1} = zeros(number_of_particles)
	A13::Array{Float64, 1} = zeros(number_of_particles)
	A21::Array{Float64, 1} = zeros(number_of_particles)
	A22::Array{Float64, 1} = zeros(number_of_particles)
	A23::Array{Float64, 1} = zeros(number_of_particles)
	A31::Array{Float64, 1} = zeros(number_of_particles)
	A32::Array{Float64, 1} = zeros(number_of_particles)
	A33::Array{Float64, 1} = zeros(number_of_particles)

	a11::Float64 = 0.0
	a12::Float64 = 0.0
	a13::Float64 = 0.0
	a21::Float64 = 0.0
	a22::Float64 = 0.0
	a23::Float64 = 0.0
	a31::Float64 = 0.0
	a32::Float64 = 0.0
	a33::Float64 = 0.0
	if particle_type != "sphere"
		for current_particle = 1:number_of_particles
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
			A11[current_particle] = a11
			A12[current_particle] = a12
			A13[current_particle] = a13
			A21[current_particle] = a21
			A22[current_particle] = a22
			A23[current_particle] = a23
			A31[current_particle] = a31
			A32[current_particle] = a32
			A33[current_particle] = a33
		end
	end

	# Create cell lists.
	cell_overlap::Float64 = 0.0
	cell_lists::Array{Array{Int64, 1}, 3} = generate_cell_lists(particle_type,
																R,
																Lx,
																Ly,
																Lz,
																X,
																Y,
																Z,
																Q0,
																Q1,
																Q2,
																Q3,
																number_of_cells_x,
																number_of_cells_y,
																number_of_cells_z,
																cell_overlap)

	# Run characterization.
	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.
	number_of_samples_per_worker::Array{Int64, 1} = convert(Array{Int64, 1}, floor(number_of_samples / number_of_workers) * ones(number_of_workers))
	number_of_samples_remaining::Int64 = number_of_samples - sum(number_of_samples_per_worker)
	number_of_samples_per_worker[1:number_of_samples_remaining] += 1
	#println(number_of_samples_per_worker)

	# Compute S2.
	S2::Array{Float64, 1} = zeros(size(d2))
#	S2 = @parallel (+) for current_worker = 1:number_of_workers
#		characterize_S2(
#			particle_type,
#			R,
#			Lx,
#			Ly,
#			Lz,
#			X,
#			Y,
#			Z,
#			Q0,
#			Q1,
#			Q2,
#			Q3,
#			A11,
#			A12,
#			A13,
#			A21,
#			A22,
#			A23,
#			A31,
#			A32,
#			A33,
#			number_of_samples_per_worker[current_worker],
#			d2,
#			cell_lists)
#	end
#	S2 /= convert(Float64, number_of_samples)

	# Compute S3.
	number_of_distances::Int64 = length(d3)
	number_of_angles::Int64 = length(theta3)
	S3::Array{Float64, 3} = zeros(number_of_distances, number_of_distances, convert(Int64, number_of_angles/2 + 1))
	S3 = @parallel (+) for current_worker = 1:number_of_workers
		characterize_S3( particle_type, R, Lx, Ly, Lz, X, Y, Z, Q0, Q1, Q2, Q3,
			A11, A12, A13, A21, A22, A23, A31, A32, A33,
			number_of_samples_per_worker[current_worker], d3, theta3, cell_lists)
	end
	S3 /= convert(Float64, number_of_samples * number_of_angles)

	# Compute zeta2.
	S1::Float64 = 1.0 - phi
	integrand::Array{Float64, 3} = zeros(number_of_distances, number_of_distances, convert(Int64, number_of_angles/2 + 1))
    for current_r1 = 1:number_of_distances
        for current_r2 = 1:number_of_distances
            for current_angle = 1:convert(Int64, number_of_angles/2 + 1)
                integrand[current_r1, current_r2, current_angle] = 1.0 / d3[current_r1] * 1.0 / d3[current_r2] *
					0.5 * sin(theta3[current_angle]) * (3.0 * cos(theta3[current_angle])^2 - 1.0) *
                    ( S3[current_r1, current_r2, current_angle] - S3[current_r1, current_r1, 1] * S3[current_r2, current_r2, 1] / S1 )
            end
        end
    end
	integral::Float64 = 0.0
	dr1::Float64 = 0.0
	dr2::Float64 = 0.0
	dtheta::Float64 = 0.0
	for current_r1 = 1:number_of_distances
		if current_r1 == 1
			dr1 = 0.5 * (d3[2] - d3[1])
		elseif current_r1 == number_of_distances
			dr1 = 0.5 * (d3[number_of_distances] - d3[number_of_distances - 1])
		else
			dr1 = 0.5 * (d3[current_r1 + 1] - d3[current_r1 - 1])
		end

        for current_r2 = 1:number_of_distances
			if current_r2 == 1
				dr2 = 0.5 * (d3[2] - d3[1])
			elseif current_r2 == number_of_distances
				dr2 = 0.5 * (d3[number_of_distances] - d3[number_of_distances - 1])
			else
				dr2 = 0.5 * (d3[current_r2 + 1] - d3[current_r2 - 1])
			end

            for current_angle = 1:convert(Int64, number_of_angles/2 + 1)
				if current_angle == 1
					dtheta = 0.5 * (theta3[2] - theta3[1])
				elseif current_angle == number_of_distances
					dtheta = 0.5 * (theta3[convert(Int64, number_of_angles/2 + 1)] - theta3[convert(Int64, number_of_angles/2 + 1) - 1])
				else
					dtheta = 0.5 * (theta3[current_angle + 1] - theta3[current_angle - 1])
				end

				integral += integrand[current_r1, current_r2, current_angle] * dr1 * dr2 * dtheta
			end
		end
	end
    zeta2::Float64 = 1.0 - 9.0/(2.0 * phi * (1.0 - phi)) * integral
	println(zeta2)

	# Kill all workers.
	#rmprocs(workers(); waitfor = typemax(Int))

	t_finish_ns::Int64 = convert(Int64, time_ns())
	t_exec::Float64 = convert(Float64, t_finish_ns - t_start_ns) / 1e9

	# Write output.
	write_output(
		output_file_path,
		number_of_samples,
		d2,
		S2,
		d3,
		theta3,
		S3,
		zeta2,
		t_exec)
	println(join(("Output written to ", output_file_path, ".")))
	println("Finished.")

	nothing
end

run_characterization()
