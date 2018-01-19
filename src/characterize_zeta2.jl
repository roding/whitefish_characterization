function characterize_zeta2(particle_type::String,
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
							cell_lists::Array{Array{Int64, 1}, 3})

	# Number of cells.
	(number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64) = size(cell_lists)

	current_particle_in_cell::Int64 = 0
	is_void::Bool = true

	#x0::Float64 = 0.0
	#y0::Float64 = 0.0
	#z0::Float64 = 0.0

	number_of_points::Int64 = 7
	x::Array{Float64, 1} = zeros(number_of_points)
	y::Array{Float64, 1} = zeros(number_of_points)
	z::Array{Float64, 1} = zeros(number_of_points)


	r::Float64 = 0.0
	s::Float64 = 0.0
	theta::Float64 = 0.0
	mu::Float64 = 1.0

	S2r::Float64 = 0.0
	S2s::Float64 = 0.0
	S3::Float64 = 0.0
	indicator::Array{Bool, 1} = falses(number_of_points)

	zeta2_integral::Array{Float64, 1} = zeros(2)
	integrand::Float64 = 0.0

	vxr::Float64 = 0.0
	vyr::Float64 = 0.0
	vzr::Float64 = 0.0
	vxs::Float64 = 0.0
	vys::Float64 = 0.0
	vzs::Float64 = 0.0
	magnitude::Float64 = 0.0

	current_cell_x::Int64 = 0
	current_cell_y::Int64 = 0
	current_cell_z::Int64 = 0
	number_of_particles_current_cell::Int64 = 0

	for current_sample = 1:number_of_samples
		#if mod(current_sample, 1000) == 0
		#	println(current_sample)
		#end

		# Pick random intial sample position for S2r.
		x[1] = Lx * rand()
		y[1] = Ly * rand()
		z[1] = Lz * rand()

		# Pick random intial sample position for S2s.
		x[3] = Lx * rand()
		y[3] = Ly * rand()
		z[3] = Lz * rand()

		# Pick random intial sample position for S3.
		x[5] = Lx * rand()
		y[5] = Ly * rand()
		z[5] = Lz * rand()

		# Pick random radii.
		r = mu * randexp()
		s = mu * randexp()

		# Pick random direction and position for r in S2r.
		vxr = randn()
		vyr = randn()
		vzr = randn()
		magnitude = sqrt(vxr * vxr + vyr * vyr + vzr * vzr)
		vxr /= magnitude
		vyr /= magnitude
		vzr /= magnitude
		x[2] = x[1] + r * vxr
		y[2] = y[1] + r * vyr
		z[2] = z[1] + r * vzr

		# Pick random direction and position for s in S2s.
		vxs = randn()
		vys = randn()
		vzs = randn()
		magnitude = sqrt(vxs * vxs + vys * vys + vzs * vzs)
		vxs /= magnitude
		vys /= magnitude
		vzs /= magnitude
		x[4] = x[3] + s * vxs
		y[4] = y[3] + s * vys
		z[4] = z[3] + s * vzs

		# Pick random 'r' sample position in random direction for S3.
		vxr = randn()
		vyr = randn()
		vzr = randn()
		magnitude = sqrt(vxr * vxr + vyr * vyr + vzr * vzr)
		vxr /= magnitude
		vyr /= magnitude
		vzr /= magnitude
		x[6] = x[5] + r * vxr
		y[6] = y[5] + r * vyr
		z[6] = z[5] + r * vzr

		# Pick random 's' sample position in random direction for S3.
		vxs = randn()
		vys = randn()
		vzs = randn()
		magnitude = sqrt(vxs * vxs + vys * vys + vzs * vzs)
		vxs /= magnitude
		vys /= magnitude
		vzs /= magnitude
		x[7] = x[5] + s * vxs
		y[7] = y[5] + s * vys
		z[7] = z[5] + s * vzs

		# Compute angle between the two vectors.
		theta = acos(vxs * vxr + vys * vyr + vzs * vzr)

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

		# Compute S2 & S3.
		S2r = 0.0
		S2s = 0.0
		S3 = 0.0

		#integrand = 0.0

		if indicator[1] & indicator[2]
			S2r = 1.0
		end
		if indicator[3] & indicator[4]
			S2s = 1.0
		end
		if indicator[5] & indicator[6] & indicator[7]
			S3 = 1.0
		end

		integrand1 = 0.5 * (3.0 * cos(theta)^2 - 1.0) * S3
		integrand2 = 0.5 * (3.0 * cos(theta)^2 - 1.0) * (- S2r * S2s)

		integrand1 /= 1.0 / mu^2 * exp(-r/mu) * exp(-s/mu) * 0.5 # This is the 'weighting' compensation in importance sampling.
		integrand2 /= 1.0 / mu^2 * exp(-r/mu) * exp(-s/mu) * 0.5 # This is the 'weighting' compensation in importance sampling.

		zeta2_integral[1] += integrand1
		zeta2_integral[2] += integrand2

	end

	return zeta2_integral
end
