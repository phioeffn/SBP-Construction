mutable struct FSBP
	xL		# Left and right ends
	xR
	N
	cpoints		# N array of collocation points
	bwidth		# allowed bandwith of the derivative matrix
	funcs		# list of basis functions
	dxfuncs		# derivatives of basis functions
	dob		# discretely orthogonal basis
	dobx		# discretely orthogonal basis derivative
	Dx
	P
	Pi
	Bx
end

mutable struct LinAdv
	xL		# Left and right domain edges
	xR
	Grid
	dx		# Size of the reference element in physical space
	dt		# timestep size
	Op		# used operator
	ul		# used BCs
	ur		#
	N
end
