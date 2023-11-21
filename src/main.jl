

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
    hx = 0.1
    hy = 0.1
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


function solve(hx, hy, tau, Re; maxiter=10)

    @show Nx, Ny = grid_num(hx, hy, tau)

    y = zeros(2 * (Nx + 1) * (Ny + 1))
    x = deepcopy(y)

    A, b = matrix(x, y, hx, hy, Nx, Ny, tau, Re)

    println("Matrix A -> size: $(size(A))")
    println("Matrix b -> size: $(size(b))")

    error = 1
    iter = 1
    error_save = []
    while error > 1e-6 && iter <= maxiter
        # if disp_error
        #     print('iter: {}'.format(iter))
        # end

        if iter == 1
            A, b = matrix(x, x, hx, hy, Nx, Ny, tau, Re)
            x = A \ b
            error = norm(x - y)
            @info "iter: " * lpad("$iter ", 3) * "error: $(error)"
            iter += 1
        else
            A, b = matrix(x, y, hx, hy, Nx, Ny, tau, Re)

            # save time step n to y
            y = x

            # solve linear sys
            x = A \ b

            # calculate error
            u = x[2:2:end]
            v = x[1:2:end]
            error = maximum(abs.((x-y)[2:end-1]))

            @info "iter: " * lpad("$iter ", 3) * "error: $(error)"


            append!(error_save, error)

            # if disp_error
            #     print('max   abs u: {}'.format(np.amax(np.abs(u))))
            #     print('max   abs v: {}'.format(np.amax(np.abs(v))))
            #     print('error abs v: {}'.format(error))
            # end
            iter += 1
        end
    end

    # if disp_error == false
    #     print('iter: {}'.format(iter))
    # end

    # # elapsed time
    # elapsed_time = time.time() - start_time
    # print('time = ', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    # # u is even elements of x
    # u = [x[i] for i in range(0, 2 * (Nx + 1) * (Ny + 1), 2)]
    u = reshape(u, (Nx + 1, Ny + 1))

    # # v is odd elements of x
    # v = [x[i] for i in range(1, 2 * (Nx + 1) * (Ny + 1), 2)]
    v = reshape(v, (Nx + 1, Ny + 1))

    # # this for export to CSV
    # if save
    #     (X, Y) = grid(hx, hy, r=1)
    #     np.savetxt("X_stream_{}_{}_{}.csv".format(Nx, Ny, Re), X, delimiter=",")
    #     np.savetxt("Y_stream_{}_{}_{}.csv".format(Nx, Ny, Re), Y, delimiter=",")
    #     np.savetxt("xi_{}_{}_{}.csv".format(Nx, Ny, Re), u, delimiter=",")
    #     np.savetxt("phi_{}_{}_{}.csv".format(Nx, Ny, Re), v, delimiter=",")
    # end

    elapsed_time = 10
    # error_save = 10

    return u, v, elapsed_time, error_save, iter
end