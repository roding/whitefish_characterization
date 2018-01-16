function characterize(	particle_type::String,
					R::Array{Float64, 2},
					Lx::Float64,
					Ly::Float64,
					Lz::Float64,
					X::Array{Float64, 1},
					Y::Array{Float64, 1},
					Z::Array{Float64, 1},
					Q0::Array{Float64, 1},
					Q1::Array{Float64, 1},
					Q2::Array{Float64, 1},
					Q3::Array{Float64, 1},
					A11::Array{Float64, 1},
					A12::Array{Float64, 1},
					A13::Array{Float64, 1},
					A21::Array{Float64, 1},
					A22::Array{Float64, 1},
					A23::Array{Float64, 1},
					A31::Array{Float64, 1},
					A32::Array{Float64, 1},
					A33::Array{Float64, 1},
					number_of_samples::Int64,
					d::Array{Float64, 1},
					cell_lists::Array{Array{Int64, 1}, 3})

	# Number of cells.
	(number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64) = size(cell_lists)

	current_particle_in_cell::Int64 = 0
	is_void::Bool = true

	number_of_points::Int64 = length(d)
	x::Array{Float64, 1} = zeros(number_of_points)
	y::Array{Float64, 1} = zeros(number_of_points)
	z::Array{Float64, 1} = zeros(number_of_points)
	S2::Array{Float64, 1} = zeros(number_of_points)
	indicator::Array{Bool, 1} = falses(number_of_points)

	vx::Float64 = 0.0
	vy::Float64 = 0.0
	vz::Float64 = 0.0

	current_cell_x::Int64 = 0
	current_cell_y::Int64 = 0
	current_cell_z::Int64 = 0
	number_of_particles_current_cell::Int64 = 0

	for current_sample = 1:number_of_samples
		if mod(current_sample, 100) == 0
			println(current_sample)
		end

		# Pick random intial sample position.
		x0 = Lx * rand()
		y0 = Ly * rand()
		z0 = Lz * rand()

		# Pick random direction vector and normalize.
		vx = randn()
		vy = randn()
		vz = randn()
		magnitude = sqrt(vx * vx + vy * vy + vz * vz)
		vx /= magnitude
		vy /= magnitude
		vz /= magnitude

		# Create vectors of sample coordinates in relative coordinates.
		x = x0 + vx * d
		y = y0 + vy * d
		z = z0 + vz * d

		# Convert to absolute coordinates (w.r.t. to Lx, Ly, Lz)
		for current_point = 1:number_of_points
			x[current_point] = position_mod(x[current_point], Lx)
			y[current_point] = position_mod(y[current_point], Ly)
			z[current_point] = position_mod(z[current_point], Lz)
		end

		# Check whether void or solid for each sample point.
		for current_point = 1:number_of_points

			current_cell_x = convert(Int64, ceil(x[current_point] / Lx * convert(Float64, number_of_cells_x)))
			current_cell_y = convert(Int64, ceil(y[current_point] / Ly * convert(Float64, number_of_cells_y)))
			current_cell_z = convert(Int64, ceil(z[current_point] / Lz * convert(Float64, number_of_cells_z)))
			number_of_particles_current_cell = length(cell_lists[current_cell_x, current_cell_y, current_cell_z])
			current_cell_list = cell_lists[current_cell_x, current_cell_y, current_cell_z]

			is_void = true
			current_particle_in_cell = 0
			while current_particle_in_cell < number_of_particles_current_cell && is_void
				current_particle_in_cell += 1

				vx = signed_distance_mod(x[current_point], X[current_cell_list[current_particle_in_cell]], Lx)
				vy = signed_distance_mod(y[current_point], Y[current_cell_list[current_particle_in_cell]], Ly)
				vz = signed_distance_mod(z[current_point], Z[current_cell_list[current_particle_in_cell]], Lz)

				if particle_type == "sphere"
					if vx^2 + vy^2 + vz^2 <= R[current_cell_list[current_particle_in_cell], 1]^2
						is_void = false
					end
				elseif particle_type == "ellipse"
					# Not supported for S2.
				elseif particle_type == "ellipsoid"
					if vx * (A11[current_cell_list[current_particle_in_cell]] * vx + A12[current_cell_list[current_particle_in_cell]] * vy + A13[current_cell_list[current_particle_in_cell]] * vz) + vy * (A21[current_cell_list[current_particle_in_cell]] * vx + A22[current_cell_list[current_particle_in_cell]] * vy + A23[current_cell_list[current_particle_in_cell]] * vz) + vz * (A31[current_cell_list[current_particle_in_cell]] * vx + A32[current_cell_list[current_particle_in_cell]] * vy + A33[current_cell_list[current_particle_in_cell]] * vz) <= 1.0
						is_void = false
					end
				elseif particle_type == "cuboid"
					(vx, vy, vz) = (A11[current_cell_list[current_particle_in_cell]] * vx + A12[current_cell_list[current_particle_in_cell]] * vy + A13[current_cell_list[current_particle_in_cell]] * vz,
									A21[current_cell_list[current_particle_in_cell]] * vx + A22[current_cell_list[current_particle_in_cell]] * vy + A23[current_cell_list[current_particle_in_cell]] * vz,
									A31[current_cell_list[current_particle_in_cell]] * vx + A32[current_cell_list[current_particle_in_cell]] * vy + A33[current_cell_list[current_particle_in_cell]] * vz)

					if abs(vx) <= R[current_cell_list[current_particle_in_cell], 1] && abs(vy) <= R[current_cell_list[current_particle_in_cell], 2] && abs(vz) <= R[current_cell_list[current_particle_in_cell], 3]
						is_void = false
					end
				elseif particle_type == "superellipsoid"
					(vx, vy, vz) = (A11[current_cell_list[current_particle_in_cell]] * vx + A12[current_cell_list[current_particle_in_cell]] * vy + A13[current_cell_list[current_particle_in_cell]] * vz,
									A21[current_cell_list[current_particle_in_cell]] * vx + A22[current_cell_list[current_particle_in_cell]] * vy + A23[current_cell_list[current_particle_in_cell]] * vz,
									A31[current_cell_list[current_particle_in_cell]] * vx + A32[current_cell_list[current_particle_in_cell]] * vy + A33[current_cell_list[current_particle_in_cell]] * vz)

					if (abs(vx)/R[current_cell_list[current_particle_in_cell], 1])^R[current_cell_list[current_particle_in_cell], 4] +
						(abs(vy)/R[current_cell_list[current_particle_in_cell], 2])^R[current_cell_list[current_particle_in_cell], 4] +
						(abs(vz)/R[current_cell_list[current_particle_in_cell], 3])^R[current_cell_list[current_particle_in_cell], 4] <= 1.0
						is_void = false
					end
				end
			end

			if is_void
				indicator[current_point] = true
			else
				indicator[current_point] = false
			end
		end

		if indicator[1] # If the origin point is in void, do something, otherwise don't.
			for current_point = 1:number_of_points
				if indicator[current_point]
					S2[current_point] += 1.0
				end
			end
		end
	end

	return S2
end
