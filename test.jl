using SparseArrays, LinearAlgebra, Plots

xr = 0
xl = 1

yr = 0
yl = 1

hx = 1 / 100
hy = 1 / 100

# the number
Nx = Int64((xl - xr) / (hx))
Ny = Int64((yl - yr) / (hy))

# function return the position of u and in vector
function index_u(i, j)
    if i > Nx + 1 || j > Ny + 1
        return false
    else
        return 2 * ((j - 1) * (Ny + 1) + i) - 1
    end
end

function index_v(i, j)
    if i > Nx + 1 || j > Ny + 1
        return false
    else
        return 2 * ((j - 1) * (Ny + 1) + i)
    end
end
# create empty vector for storing element in matrix
index_i = []
index_j = []

function sys_stream(x::Array, y::Array; tau=hx, Re=0)
    index_i = []
    index_j = []
    index_k = []

    # create empty right hamd side vector
    b = zeros(2 * (Nx + 1) * (Ny + 1))

    epsilon = 1e-8
    for i in 1:Nx + 1
        for j in 1:Ny + 1
            if i == Nx + 1 && j == Ny + 1

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, -1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, hy * epsilon)

                b[index_u(i, j)] = hy
                b[index_v(i, j)] = 0

            elseif i == 1 && j == Ny + 1

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, -1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, hy * epsilon)

                b[index_u(i, j)] = hy
                b[index_v(i, j)] = 0

            elseif i == 1 && j == 1

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -2)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, hy * epsilon)

                b[index_u(i, j)] = hy
                b[index_v(i, j)] = 0

            elseif j == Ny + 1
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, -1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, hy * epsilon)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                b[index_u(i, j)] = hy
                b[index_v(i, j)] = 0
            elseif i == 1

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, epsilon)

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0
            elseif i == Nx + 1
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i - 1, j))
                append!(index_k, -1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, epsilon)

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0

            elseif j == 1
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, epsilon)

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0
            else
                # nonlinear terms
                phi_x = x[index_v(i + 1, j)] - x[index_v(i - 1, j)]
                phi_y = x[index_v(i, j + 1)] - x[index_v(i, j - 1)]

                xi_x = x[index_u(i + 1, j)] - x[index_u(i - 1, j)]
                xi_y = x[index_u(i, j + 1)] - x[index_u(i, j - 1)]

                # right hand side
                # first equation
                b[index_u(i, j)] = x[index_u(i, j)] + tau / (2 * hx^2) * (
                    x[index_u(i + 1, j)] - 4 * x[index_u(i, j)] +
                    x[index_u(i - 1, j)] + x[index_u(i, j + 1)] +
                    x[index_u(i, j - 1)])

                # second equation
                b[index_v(i, j)] = -(x[index_u(i, j)])

                # diagonal of matrix
                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, 1 + 2. * tau / (hx^2))

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -4 / (hx^2))

                
                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j + 1))
                append!(index_k, tau * Re / (8 * hx * hy) * phi_x - tau /
                    (2. * hx^2))

                
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, -tau * Re / (8 * hx * hy) * xi_x)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j - 1))
                append!(index_k, -tau * Re / (8 * hx * hy) * phi_x - tau /
                    (2. * hx^2))

                
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, tau * Re / (8 * hx * hy) * xi_x)

                
                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i + 1, j))
                append!(index_k, -tau * Re / (8 * hx * hy) * phi_y - tau /
                    (2. * hx^2))

                
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, tau * Re / (8 * hx * hy) * xi_y)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i - 1, j))
                append!(index_k, tau * Re / (8 * hx * hy) * phi_y - tau /
                    (2. * hx^2))

                
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i - 1, j))
                append!(index_k, -tau * Re / (8 * hx * hy) * xi_y)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, 1 / (hy^2))

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, 1 / (hy^2))

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, 1 / (hx^2))

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i - 1, j))
                append!(index_k, 1 / (hx^2))

                append!(index_i, index_v(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, 1)
            end
        end
    end
    return (sparse(index_i, index_j, index_k), b)
end

y = zeros(2 * (Nx + 1) * (Ny + 1))

# Loop updating x from x=A^-1*b
# x is solution at time step n+1
# y is solution at time step n
x = copy(y)
error = 1
iteration = 1
error_save = []
max_iteration = 150
tau = hx
Re = 1000

# @gif for i in 1:max_iteration
while error > 1e-8 && iteration <= max_iteration
    global iteration, x, y, error_save, u, v, error
    if iteration == 1
        A, b = sys_stream(x, x, tau=tau, Re=Re)
        x = A \ b
        error = norm(x - y)
        u = [x[i] for i in 1:2:2 * (Nx + 1) * (Ny + 1)]
        v = [x[i] for i in 2:2:2 * (Nx + 1) * (Ny + 1)]
        iteration += 1
    else
        A, b = sys_stream(x, y, tau=tau, Re=Re)
        
        # save time step n to y
        y = copy(x)
        
        # solve linear sys
        x = A \ b

        # calculate error
        u = [x[i] for i in 1:2:2 * (Nx + 1) * (Ny + 1)]
        v = [x[i] for i in 2:2:2 * (Nx + 1) * (Ny + 1)]
        error = maximum(abs.([x[i] - y[i] for i in 2:2:2 * (Nx + 1) * (Ny + 1)]))
        append!(error_save, error)
        println("iteration: $(iteration)")
        println("error    : $(error)")
        iteration += 1
    end
end


u = [x[i] for i in 1:2:2 * (Nx + 1) * (Ny + 1)]
v = [x[i] for i in 2:2:2 * (Nx + 1) * (Ny + 1)]
u = reshape(u, (Nx + 1, Ny + 1))
v = reshape(v, (Nx + 1, Ny + 1))

phi = v
value(x, y) = phi[x, y]
m = minimum(phi)
M = maximum(phi)
sc1 = [i for i in m:-m/10:0]
sc2 = [i for i in 0:M/10:M]
sc = append!(sc1, sc2)
Plots.contour(1:Nx, 1:Ny, value, levels=sc, c=:leonardo)
# Plots.contour(1:52, 1:52, value, fill=true)
# m = minimum(phi)
# M = maximum(phi)
# sc1 = [i for i in m:-m/10:0]
# sc2 = [i for i in 0:M/10:M]
# sc = append!(sc1, sc2)
# Plots.contour(1:52, 1:52, value, levels=sc)
# sc = [m:-m/20:0, 0:M/20:M]
# cp = plt.contour(phi, levels=sc, cmap='RdGy')
# using GR
# GR.contourf(phi, levels=sc)
# Plots.contour(1:52, 1:52, value, fill=true) 
# Plots.heatmap(1:52, 1:52, value)
# GR.heatmap(phi)