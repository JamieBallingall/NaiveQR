using LinearAlgebra

# Compute the eigenvalues of a matrix using the QR algorithm
function eigvals_qr(A :: Matrix{Float64}, tol :: Float64 = 1e-9, maxiter :: Int64 = 10000) :: Vector{Float64}
  i = 0
  while i < maxiter && worst_off_diag(A) > tol
    (Q, R) = qr(A)
    A = R * Q
    i += 1
  end
  return i >= maxiter ? Vector{Float64}() : sort(diag(A))
end

# Return the absolute value of the largest (in absolute value) off diagonal element of a matrix
function worst_off_diag(A :: Matrix{Float64}) :: Float64
  return maximum(abs.(A .- diagm(size(A, 1), size(A, 2), diag(A))))
end
