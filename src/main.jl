

function grid_num(hx::Real, hy::Real, tau::Real)
    # boundrayy
    xr, xl = 0, 1
    yr, yl = 0, 1

    # calculate number of create_grid
    Nx = Int64((xl-xr)/hx)
    Ny = Int64((yl-yr)/hy)

    return Nx, Ny

end


function matrix()
    hx = 0.001
    hy = 0.001
    tau = 0.1
    Re = 100
    Nx, Ny = grid_num(hx, hy, tau)
    y = zeros(2 * (Nx + 1) * (Ny + 1))
    x = deepcopy(y)

    A, ~ = matrix(x, y, hx, hy, Nx, Ny, tau, Re)
    return A
end


function matrix(x, y, hx, hy, Nx, Ny, tau, Re)

    # parameter
    epsilon1 = 1e-8  # 1 =  Nx, Ny
    epsilon2 = 1e-8  # 2 =   0, Ny
    epsilon3 = 1e-8  # 3 =   0, Ny
    epsilon4 = 1e-8  # 4 =   0,  0
    epsilon5 = 1e-8  # 5 =   i = Nx
    epsilon6 = 1e-8  # 6 =   i = 0

    index_i = []
    index_j = []
    index_k = []

    # create empty right hamd side vector
    b = zeros(2 * (Nx + 1) * (Ny + 1))

    function index_u(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 1
    end
    
    function index_v(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 2
    end

    for i in 0:(Nx)
        for j in 0:(Ny)
            if i == Nx && j == Ny

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i - 1, j))
                append!(index_k, -1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                # perturb
                if epsilon1 != 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, epsilon1)
                end

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0

            elseif  i == 0 && j == Ny
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                # index_i.append!(index_v(i, j))
                # index_j.append!(index_v(i, j-1))
                # index_k.append!(1)

                if epsilon2 != 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, epsilon2)
                end

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0
            
            elseif  i == 0 && j == 0

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, 1)

                # d/dx(phi)
                # index_i.append!(index_u(i, j))
                # index_j.append!(index_v(i+1, j))
                # index_k.append!(1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, 1)

                if epsilon3 != 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, epsilon3)
                end

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0

            elseif  j == Ny

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, -1)

                # index_i.append!(index_u(i, j))
                # index_j.append!(index_u(i, j))
                # index_k.append!(epsilon)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, 1)

                b[index_u(i, j)] = hy
                b[index_v(i, j)] = 0

            elseif  i == 0

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, 1)

                if epsilon4 != 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, epsilon4)
                end

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0

            elseif  i == Nx

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i - 1, j))
                append!(index_k, -1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, 1)

                if epsilon5 != 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, epsilon5)
                end

                b[index_u(i, j)] = 0
                b[index_v(i, j)] = 0

            elseif  j == 0

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -1)

                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, 1)

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, 1)

                if epsilon6 != 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, epsilon6)
                end

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
                b[index_u(i, j)] = x[index_u(i, j)] + tau / (2 * hx^2) * ( x[index_u(i + 1, j)] - 4 * x[index_u(i, j)] + x[index_u(i - 1, j)] + x[index_u(i, j + 1)] + x[index_u(i, j - 1)])

                # second equation
                b[index_v(i, j)] = -(x[index_u(i, j)])
                b[index_v(i, j)] -= (x[index_v(i + 1, j)]) / (hx^2)
                b[index_v(i, j)] -= (x[index_v(i - 1, j)]) / (hx^2)
                b[index_v(i, j)] -= (x[index_v(i, j + 1)]) / (hx^2)
                b[index_v(i, j)] -= (x[index_v(i, j - 1)]) / (hx^2)
                b[index_v(i, j)] += 4 * (x[index_v(i, j)]) / (hx^2)

                # diagonal of matrix
                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j))
                append!(index_k, 1 + 2. * tau / (hx^2))

                append!(index_i, index_v(i, j))
                append!(index_j, index_v(i, j))
                append!(index_k, -4 / (hx^2))

                # u x
                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j + 1))
                append!(index_k, tau * Re / (8 * hx * hy) * phi_x - tau / (2. * hx^2))

                # add
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j + 1))
                append!(index_k, -tau * Re / (8 * hx * hy) * xi_x)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i, j - 1))
                append!(index_k, -tau * Re / (8 * hx * hy) * phi_x - tau / (2. * hx^2))

                # add
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i, j - 1))
                append!(index_k, tau * Re / (8 * hx * hy) * xi_x)

                # u y
                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i + 1, j))
                append!(index_k, -tau * Re / (8 * hx * hy) * phi_y - tau / (2. * hx^2))

                # add
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i + 1, j))
                append!(index_k, tau * Re / (8 * hx * hy) * xi_y)

                append!(index_i, index_u(i, j))
                append!(index_j, index_u(i - 1, j))
                append!(index_k, tau * Re / (8 * hx * hy) * phi_y - tau / (2. * hx^2))

                # add
                append!(index_i, index_u(i, j))
                append!(index_j, index_v(i - 1, j))
                append!(index_k, -tau * Re / (8 * hx * hy) * xi_y)

                # equation v
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

    return sparse(index_i, index_j, index_k, 2 * (Nx + 1) * (Ny + 1), 2 * (Nx + 1) * (Ny + 1)), b
end


function print_Error(iter, maxiter, Error)
    max_pad = length(string(maxiter))
    print("|$(Dates.format(Dates.now(), "HH:MM:SS"))|" * lpad("$iter", max_pad) * "/$maxiter|")
    @printf("%.5e|\n", Error)
end


function solve(hx, hy, tau, Re; maxiter=10, figsave=true)

    @show Nx, Ny = grid_num(hx, hy, tau)

    y = zeros(2 * (Nx + 1) * (Ny + 1))
    x = y

    A, b = matrix(x, y, hx, hy, Nx, Ny, tau, Re)

    println("Matrix A -> size: $(size(A))")
    println("Matrix b -> size: $(size(b))")

    Error = 1
    iter = 1
    Error_save = []
    while Error > 1e-6 && iter <= maxiter

        if iter == 1
            A, b = matrix(x, x, hx, hy, Nx, Ny, tau, Re)
            x = A \ b
            Error = norm(x - y)
            print_Error(iter, maxiter, Error)
            iter += 1
        else
            A, b = matrix(x, y, hx, hy, Nx, Ny, tau, Re)
            
            # save time step n to y
            y = x
            
            # solve linear sys
            # x = A \ b
            luA = lu(A)
            x = luA \ b
            
            # calculate Error
            u = x[2:2:end]
            v = x[1:2:end]
            Error = maximum(abs.((x-y)[2:end-1]))
            
            print_Error(iter, maxiter, Error)
            append!(Error_save, Error)
            iter += 1
        end
    end

    u = reshape(u, (Nx + 1, Ny + 1))
    v = reshape(v, (Nx + 1, Ny + 1))

    elapsed_time = 10

    if figsave
        # fig_name = "$(Dates.format(Dates.now(), "YY-MM-DD"))"
        fig_name = "hx_$(hx)_hy_$(hy)_tau_$(tau)_Re_$(Re)"
        savefig(contour(u, fill=true, title=fig_name), "fig/U/$(fig_name)_U.png");
        savefig(contour(v, fill=true, title=fig_name), "fig/V/$(fig_name)_V.png");
        @info "Done!!!"
    end

    # return u, v, elapsed_time, Error_save, iter
end