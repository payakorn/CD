using CD

# 
hx = 0.1
hy = hx
tau = hx
Re = 100
t = 1
y = zeros(2 * (Nx + 1) * (Ny + 1))
x = y

# 
Nx = 10
Ny = 10
print("hx = $hx, tau = $tau, t = $t")

# function
CD.generate_matrix(x, y, hx, hy, Nx, Ny, tau)
