# Minimal FFI smoke test: call Fortran `square_` from a shared library.
# Update the absolute library path below if your local checkout differs.
x = 5.0
result = Ref{Float64}(0.0)
ccall((:square_, "/calculate/Gadwala/Github_projects/ferrite-fortran-integration_using_Julia/test/libsquare.so"), Cvoid, (Ref{Float64}, Ref{Float64}), x, result)
println("Square of $x is $(result[])")