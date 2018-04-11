function characterize_S3(	particle_type::String,
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
							d3::Array{Float64, 1},
							theta3::Array{Float64, 1},
							cell_lists::Array{Array{Int64, 1}, 3})

	# Number of cells.
	(number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64) = size(cell_lists)

	current_particle_in_cell::Int64 = 0
	is_void::Bool = true

	number_of_distances::Int64 = length(d3)
	number_of_angles::Int64 = length(theta3)

	x0::Float64 = 0.0
	y0::Float64 = 0.0
	z0::Float64 = 0.0

	x_rel::Array{Float64, 2} = zeros(number_of_distances, number_of_angles)
	y_rel::Array{Float64, 2} = zeros(number_of_distances, number_of_angles)
	z_rel::Array{Float64, 2} = zeros(number_of_distances, number_of_angles)
	for current_distance = 1:number_of_distances
		for current_angle = 1:number_of_angles
			x_rel[current_distance, current_angle] = d3[current_distance] * cos(theta3[current_angle])
			y_rel[current_distance, current_angle] = d3[current_distance] * sin(theta3[current_angle])
		end
	end
	x::Array{Float64, 2} = zeros(number_of_distances, number_of_angles)
	y::Array{Float64, 2} = zeros(number_of_distances, number_of_angles)
	z::Array{Float64, 2} = zeros(number_of_distances, number_of_angles)

	indicator_original::BitArray{2} = falses(number_of_distances, number_of_angles)
	indicator::BitArray{2} = falses(number_of_distances, number_of_angles)
	S3::Array{Float64, 3} = zeros(number_of_distances, number_of_distances, convert(Int64, number_of_angles/2 + 1))

	vx::Float64 = 0.0
	vy::Float64 = 0.0
	vz::Float64 = 0.0

	current_cell_x::Int64 = 0
	current_cell_y::Int64 = 0
	current_cell_z::Int64 = 0
	number_of_particles_current_cell::Int64 = 0

	q0::Float64 = 0.0
	q1::Float64 = 0.0
	q2::Float64 = 0.0
	q3::Float64 = 0.0
	a11::Float64 = 0.0
	a12::Float64 = 0.0
	a13::Float64 = 0.0
	a21::Float64 = 0.0
	a22::Float64 = 0.0
	a23::Float64 = 0.0
	a31::Float64 = 0.0
	a32::Float64 = 0.0
	a33::Float64 = 0.0

	for current_sample = 1:number_of_samples
		#if mod(current_sample, 1000) == 0
		#	println(current_sample)
		#end

		# Pick random intial sample position.
		x0 = Lx * rand()
		y0 = Ly * rand()
		z0 = Lz * rand()

		# Check whether void or solid.
		current_cell_x = convert(Int64, ceil(x0 / Lx * convert(Float64, number_of_cells_x)))
		current_cell_y = convert(Int64, ceil(y0 / Ly * convert(Float64, number_of_cells_y)))
		current_cell_z = convert(Int64, ceil(z0 / Lz * convert(Float64, number_of_cells_z)))
		number_of_particles_current_cell = length(cell_lists[current_cell_x, current_cell_y, current_cell_z])
		current_cell_list = cell_lists[current_cell_x, current_cell_y, current_cell_z]

		is_void = true
		current_particle_in_cell = 0
		while current_particle_in_cell < number_of_particles_current_cell && is_void
			current_particle_in_cell += 1

			vx = signed_distance_mod(x0, X[current_cell_list[current_particle_in_cell]], Lx)
			vy = signed_distance_mod(y0, Y[current_cell_list[current_particle_in_cell]], Ly)
			vz = signed_distance_mod(z0, Z[current_cell_list[current_particle_in_cell]], Lz)

			if particle_type == "sphere"
				if vx^2 + vy^2 + vz^2 <= R[current_cell_list[current_particle_in_cell], 1]^2
					is_void = false
				end
			else
				(vx, vy, vz) = (A11[current_cell_list[current_particle_in_cell]] * vx + A12[current_cell_list[current_particle_in_cell]] * vy + A13[current_cell_list[current_particle_in_cell]] * vz,
								A21[current_cell_list[current_particle_in_cell]] * vx + A22[current_cell_list[current_particle_in_cell]] * vy + A23[current_cell_list[current_particle_in_cell]] * vz,
								A31[current_cell_list[current_particle_in_cell]] * vx + A32[current_cell_list[current_particle_in_cell]] * vy + A33[current_cell_list[current_particle_in_cell]] * vz)
				if particle_type == "ellipsoid"
					if  (vx/R[current_cell_list[current_particle_in_cell], 1])^2 +
						(vy/R[current_cell_list[current_particle_in_cell], 2])^2 +
						(vz/R[current_cell_list[current_particle_in_cell], 3])^2 <= 1.0
						is_void = false
					end
				elseif particle_type == "cuboid"
					if abs(vx) <= R[current_cell_list[current_particle_in_cell], 1] && abs(vy) <= R[current_cell_list[current_particle_in_cell], 2] && abs(vz) <= R[current_cell_list[current_particle_in_cell], 3]
						is_void = false
					end
				elseif particle_type == "superellipsoid"
					if  (abs(vx)/R[current_cell_list[current_particle_in_cell], 1])^R[current_cell_list[current_particle_in_cell], 4] +
						(abs(vy)/R[current_cell_list[current_particle_in_cell], 2])^R[current_cell_list[current_particle_in_cell], 4] +
						(abs(vz)/R[current_cell_list[current_particle_in_cell], 3])^R[current_cell_list[current_particle_in_cell], 4] <= 1.0
						is_void = false
					end
				end
			end
		end

		# If origin position is in void, continue analysis.
		if is_void
			# Pick a random orientation.
			(q0, q1, q2, q3) = generate_random_unit_quaternion()
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(q0, q1, q2, q3)

			# Rotate the relative coordinates of the 'template' (as adapted from Hlushkou, 2015)
			# and add the origin coordinate.
			for current_distance = 1:number_of_distances
				for current_angle = 1:number_of_angles
					x[current_distance, current_angle] = x0 + a11 * x_rel[current_distance, current_angle] + a12 * y_rel[current_distance, current_angle] + a13 * z_rel[current_distance, current_angle]
					y[current_distance, current_angle] = y0 + a21 * x_rel[current_distance, current_angle] + a22 * y_rel[current_distance, current_angle] + a23 * z_rel[current_distance, current_angle]
					z[current_distance, current_angle] = z0 + a31 * x_rel[current_distance, current_angle] + a32 * y_rel[current_distance, current_angle] + a33 * z_rel[current_distance, current_angle]
				end
			end

			# Convert to absolute coordinates (w.r.t. to Lx, Ly, Lz)
			for current_distance = 1:number_of_distances
				for current_angle = 1:number_of_angles
					x[current_distance, current_angle] = position_mod(x[current_distance, current_angle], Lx)
					y[current_distance, current_angle] = position_mod(y[current_distance, current_angle], Ly)
					z[current_distance, current_angle] = position_mod(z[current_distance, current_angle], Lz)
				end
			end

			# Check whether void or solid for each sample point.
			for current_distance = 1:number_of_distances
				for current_angle = 1:number_of_angles
					current_cell_x = convert(Int64, ceil(x[current_distance, current_angle] / Lx * convert(Float64, number_of_cells_x)))
					current_cell_y = convert(Int64, ceil(y[current_distance, current_angle] / Ly * convert(Float64, number_of_cells_y)))
					current_cell_z = convert(Int64, ceil(z[current_distance, current_angle] / Lz * convert(Float64, number_of_cells_z)))
					number_of_particles_current_cell = length(cell_lists[current_cell_x, current_cell_y, current_cell_z])
					current_cell_list = cell_lists[current_cell_x, current_cell_y, current_cell_z]

					is_void = true
					current_particle_in_cell = 0
					while current_particle_in_cell < number_of_particles_current_cell && is_void
						current_particle_in_cell += 1

						vx = signed_distance_mod(x[current_distance, current_angle], X[current_cell_list[current_particle_in_cell]], Lx)
						vy = signed_distance_mod(y[current_distance, current_angle], Y[current_cell_list[current_particle_in_cell]], Ly)
						vz = signed_distance_mod(z[current_distance, current_angle], Z[current_cell_list[current_particle_in_cell]], Lz)

						if particle_type == "sphere"
							if vx^2 + vy^2 + vz^2 <= R[current_cell_list[current_particle_in_cell], 1]^2
								is_void = false
							end
						else
							(vx, vy, vz) = (A11[current_cell_list[current_particle_in_cell]] * vx + A12[current_cell_list[current_particle_in_cell]] * vy + A13[current_cell_list[current_particle_in_cell]] * vz,
											A21[current_cell_list[current_particle_in_cell]] * vx + A22[current_cell_list[current_particle_in_cell]] * vy + A23[current_cell_list[current_particle_in_cell]] * vz,
											A31[current_cell_list[current_particle_in_cell]] * vx + A32[current_cell_list[current_particle_in_cell]] * vy + A33[current_cell_list[current_particle_in_cell]] * vz)
							if particle_type == "ellipsoid"
								if  (vx/R[current_cell_list[current_particle_in_cell], 1])^2 +
									(vy/R[current_cell_list[current_particle_in_cell], 2])^2 +
									(vz/R[current_cell_list[current_particle_in_cell], 3])^2 <= 1.0
									is_void = false
								end
							elseif particle_type == "cuboid"
								if abs(vx) <= R[current_cell_list[current_particle_in_cell], 1] && abs(vy) <= R[current_cell_list[current_particle_in_cell], 2] && abs(vz) <= R[current_cell_list[current_particle_in_cell], 3]
									is_void = false
								end
							elseif particle_type == "superellipsoid"
								if  (abs(vx)/R[current_cell_list[current_particle_in_cell], 1])^R[current_cell_list[current_particle_in_cell], 4] +
									(abs(vy)/R[current_cell_list[current_particle_in_cell], 2])^R[current_cell_list[current_particle_in_cell], 4] +
									(abs(vz)/R[current_cell_list[current_particle_in_cell], 3])^R[current_cell_list[current_particle_in_cell], 4] <= 1.0
									is_void = false
								end
							end
						end
					end

					if is_void
						indicator_original[current_distance, current_angle] = true
					else
						indicator_original[current_distance, current_angle] = false
					end
				end
			end

			for current_shift = 0:number_of_angles-1
                indicator = circshift(indicator_original, (0, current_shift))
                for current_r1 = 1:number_of_distances
                    if indicator[current_r1, 1]
						for current_r2 = 1:number_of_distances
							for current_angle = 1:convert(Int64, number_of_angles / 2 + 1)
								if indicator[current_r2, current_angle]
									S3[current_r1, current_r2, current_angle] += 1.0
								end
							end
						end
                    end
                end
            end
		end
	end

	return S3
end
