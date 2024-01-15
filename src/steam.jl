# new section for Stream function
function sys_stream1(hx, hy, tau, Re; max_iteration=100, disp_error=false, epsilon=nothing, left=0, right=0, top=1, below=0, save_every=20, initial=nothing, iteration=1)

    # delete folder checkpoint
    location = "checkpoints/6_epsilon/$(hx)_$(hy)_$(tau)_$(Re)_$(epsilon[1])"
    if iteration == 1
        rm(location, recursive=true, force=true)
        mkpath(location)
    end

    # boundrary
    xr = 0
    xl = 1

    yr = 0
    yl = 1

    # the number
    Nx = Int64((xl - xr) / (hx))
    Ny = Int64((yl - yr) / (hy))

    # function return the position of u and v in vector
    function index_u(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 1
    end

    function index_v(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 2
    end

    # create empty vector for storing elements in matrix
    index_i = []
    index_j = []

    function create_matrix_stream(x, y; Nx=false, Ny=false, tau=false, Re=Re, index_i=index_i, index_j=index_j)
        index_i = []
        index_j = []
        index_k = []

        # create empty right hamd side vector
        b = zeros(2 * (Nx + 1) * (Ny + 1))

        # epsilon
        if isnothing(epsilon)
            epsilon1 = 1e-8  # 1 =  Nx, Ny
            epsilon2 = 1e-8  # 2 =   0, Ny
            epsilon3 = 1e-8  # 3 =   0,  0
            epsilon4 = 1e-8  # 4 =   i = 0
            epsilon5 = 1e-8  # 5 =   i = Nx
            epsilon6 = 1e-8  # 6 =   j = 0
        else
            (
                epsilon1,
                epsilon2,
                epsilon3,
                epsilon4,
                epsilon5,
                epsilon6,
            ) = epsilon
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

                    b[index_u(i, j)] = hx*right
                    b[index_v(i, j)] = 0

                elseif i == 0 && j == Ny
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i, j-1))
                    # append!(index_k, 1)

                    if epsilon2 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon2)
                    end

                    b[index_u(i, j)] = hx*left
                    b[index_v(i, j)] = 0

                elseif i == 0 && j == 0

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, 1)

                    # d/dx(phi)
                    # append!(index_i, index_u(i, j))
                    # append!(index_j, index_v(i+1, j))
                    # append!(index_k, 1)

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

                    b[index_u(i, j)] = hx*below
                    b[index_v(i, j)] = 0

                elseif j == Ny
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j - 1))
                    append!(index_k, -1)

                    # append!(index_i, index_u(i, j))
                    # append!(index_j, index_u(i, j))
                    # append!(index_k, epsilon)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1)

                    b[index_u(i, j)] = hy*top
                    b[index_v(i, j)] = 0

                elseif i == 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i, j + 1))
                    # append!(index_k, 1)

                    if epsilon4 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon4)
                    end

                    b[index_u(i, j)] = hx*left
                    b[index_v(i, j)] = 0

                elseif i == Nx
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i - 1, j))
                    append!(index_k, -1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i, j - 1))
                    # append!(index_k, 1)

                    if epsilon5 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon5)
                    end

                    b[index_u(i, j)] = hx*right
                    b[index_v(i, j)] = 0

                elseif j == 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i + 1, j))
                    # append!(index_k, 1)

                    if epsilon6 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon6)
                    end

                    b[index_u(i, j)] = hx*below
                    b[index_v(i, j)] = 0

                else

                    # nonlinear terms
                    phi_x = x[index_v(i + 1, j)] - x[index_v(i - 1, j)]
                    phi_y = x[index_v(i, j + 1)] - x[index_v(i, j - 1)]

                    xi_x = x[index_u(i + 1, j)] - x[index_u(i - 1, j)]
                    xi_y = x[index_u(i, j + 1)] - x[index_u(i, j - 1)]

                    # right hand side
                    # first equation
                    b[index_u(i, j)] = x[index_u(i, j)] + tau / (2 * hx ^ 2) * (
                        x[index_u(i + 1, j)]
                        - 4 * x[index_u(i, j)]
                        + x[index_u(i - 1, j)]
                        + x[index_u(i, j + 1)]
                        + x[index_u(i, j - 1)]
                    )

                    # second equation
                    b[index_v(i, j)] = -(x[index_u(i, j)])
                    b[index_v(i, j)] -= (x[index_v(i + 1, j)]) / (hx ^ 2)
                    b[index_v(i, j)] -= (x[index_v(i - 1, j)]) / (hx ^ 2)
                    b[index_v(i, j)] -= (x[index_v(i, j + 1)]) / (hx ^ 2)
                    b[index_v(i, j)] -= (x[index_v(i, j - 1)]) / (hx ^ 2)
                    b[index_v(i, j)] += 4 * (x[index_v(i, j)]) / (hx ^ 2)

                    # diagonal of matrix
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, 1 + 2.0 * tau / (hx ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -4 / (hx ^ 2))

                    # u x
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j + 1))
                    append!(index_k, tau * Re / (8 * hx * hy) * phi_x - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, -tau * Re / (8 * hx * hy) * xi_x)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j - 1))
                    append!(index_k, -tau * Re / (8 * hx * hy) * phi_x - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j - 1))
                    append!(index_k, tau * Re / (8 * hx * hy) * xi_x)

                    # u y
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i + 1, j))
                    append!(index_k, -tau * Re / (8 * hx * hy) * phi_y - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, tau * Re / (8 * hx * hy) * xi_y)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i - 1, j))
                    append!(index_k, tau * Re / (8 * hx * hy) * phi_y - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i - 1, j))
                    append!(index_k, -tau * Re / (8 * hx * hy) * xi_y)

                    # equation v
                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, 1 / (hy ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j - 1))
                    append!(index_k, 1 / (hy ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1 / (hx ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i - 1, j))
                    append!(index_k, 1 / (hx ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, 1)
                end
            end
        end

        return (sparse(index_i, index_j, index_k, 2 * (Nx + 1) * (Ny + 1), 2 * (Nx + 1) * (Ny + 1)), b)
    end

    # start
    start_time = Dates.now()

    # y is an initial solution
    if isnothing(initial)
        y = zeros(2 * (Nx + 1) * (Ny + 1))
    else
        y = initial
        iteration = iteration + 1
    end

    # Loop updating x from x=A^-1*b
    # x is solution at time step n+1
    # y is solution at time step n
    x = deepcopy(y)
    error = 1
    error_save = []
    while error > 1e-6 && iteration <= max_iteration

        if iteration == 1
            A, b = create_matrix_stream(x, x, Nx=Nx, Ny=Ny, tau=tau, Re=Re)
            luA = lu(A)
            x = luA \ b

            @info "$(Dates.format(Dates.now(), Dates.dateformat"dd u Y -- H:M:S")) iter: $iteration/$max_iteration"
            println("\t max \t abs u: $(maximum(x[1:2:end]))")
            println("\t max \t abs v: $(maximum(x[2:2:end]))")
            iteration += 1
        else
            A, b = create_matrix_stream(x, y, Nx=Nx, Ny=Ny, tau=tau, Re=Re)

            # save time step n to y
            y = deepcopy(x)

            # solve linear sys
            luA = lu(A)
            x = luA \ b

            # calculate u, v
            u = x[1:2:end]
            v = x[2:2:end]

            # calculate error
            @info "$(Dates.format(Dates.now(), Dates.dateformat"dd u Y -- H:M:S")) iter: $iteration/$max_iteration"
            println("\t max \t abs u: $(maximum(u))")
            println("\t max \t abs v: $(maximum(v))")

            error = maximum(abs.(v-y[2:2:end]))
            println("\t error \t abs v: $(error)")
            push!(error_save, error)

            # save every 10 iterations
            if iteration % save_every == 0
                # plot_contour(reshape(v, (Nx + 1, Ny + 1)), "fig/save/checkpoint_$iteration.png")
                save_checkpoint(joinpath(location, "iter_u_$iteration.txt"), reshape(u, (Nx + 1, Ny + 1)))
                save_checkpoint(joinpath(location, "iter_v_$iteration.txt"), reshape(v, (Nx + 1, Ny + 1)))
            end
            
            iteration += 1
        end
    end
    
    u = x[1:2:end]
    v = x[2:2:end]
    save_checkpoint(joinpath(location, "iter_u_$(iteration-1).txt"), reshape(u, (Nx + 1, Ny + 1)))
    save_checkpoint(joinpath(location, "iter_v_$(iteration-1).txt"), reshape(v, (Nx + 1, Ny + 1)))

    u = reshape(u, (Nx + 1, Ny + 1))
    v = reshape(v, (Nx + 1, Ny + 1))

    stop_time = Dates.now()
    elapsed_time = Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(stop_time) - Dates.DateTime(start_time)))
    @info "elapsed time = $elapsed_time"

    return u, v, elapsed_time, error_save, iteration
end


function sys_stream2(hx, hy, tau, Re; max_iteration=100, disp_error=false, epsilon=nothing, left=0, right=0, top=1, below=0, save_every=20, initial=nothing, iteration=1)

    # delete folder checkpoints
    location = "checkpoints/8_epsilon/$(hx)_$(hy)_$(tau)_$(Re)_$(epsilon[1])"
    if iteration == 1
        rm(location, recursive=true, force=true)
        mkpath(location)
    end

    # boundrary
    xr = 0
    xl = 1

    yr = 0
    yl = 1

    # the number
    Nx = Int64((xl - xr) / (hx))
    Ny = Int64((yl - yr) / (hy))

    # function return the position of u and v in vector
    function index_u(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 1
    end

    function index_v(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 2
    end

    # create empty vector for storing elements in matrix
    index_i = []
    index_j = []

    function create_matrix_stream(x, y; Nx=false, Ny=false, tau=false, Re=Re, index_i=index_i, index_j=index_j, epsilon=nothing)
        index_i = []
        index_j = []
        index_k = []

        # create empty right hamd side vector
        b = zeros(2 * (Nx + 1) * (Ny + 1))

        # epsilon
        if isnothing(epsilon)
            epsilon1 = 1e-8  # 1 =  Nx, Ny
            epsilon2 = 1e-8  # 2 =   0, Ny
            epsilon3 = 1e-8  # 3 =   0,  0
            epsilon4 = 1e-8  # 4 =   i = 0
            epsilon12 = 1e-8  # 5 =   i = Nx
            epsilon23 = 1e-8  # 6 =   j = 0
            epsilon34 = 1e-8  # 6 =   j = 0
            epsilon14 = 1e-8  # 6 =   j = 0
        else
            (
                epsilon1,
                epsilon2,
                epsilon3,
                epsilon4,
                epsilon12,
                epsilon23,
                epsilon34,
                epsilon14,
            ) = epsilon
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
                    if epsilon3 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon3)
                    end

                    b[index_u(i, j)] = hx*right
                    b[index_v(i, j)] = 0

                elseif i == Nx && j == 0
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
                    if epsilon4 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon4)
                    end

                    b[index_u(i, j)] = hx*right
                    b[index_v(i, j)] = 0

                elseif i == 0 && j == Ny
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i, j-1))
                    # append!(index_k, 1)

                    if epsilon2 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon2)
                    end

                    b[index_u(i, j)] = hx*left
                    b[index_v(i, j)] = 0

                elseif i == 0 && j == 0

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, 1)

                    # d/dx(phi)
                    # append!(index_i, index_u(i, j))
                    # append!(index_j, index_v(i+1, j))
                    # append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1)

                    if epsilon1 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon1)
                    end

                    b[index_u(i, j)] = hx*below
                    b[index_v(i, j)] = 0

                elseif j == Ny
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j - 1))
                    append!(index_k, -1)

                    # append!(index_i, index_u(i, j))
                    # append!(index_j, index_u(i, j))
                    # append!(index_k, epsilon)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1)

                    if epsilon23 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon23)
                    end

                    b[index_u(i, j)] = hy*top
                    b[index_v(i, j)] = 0

                elseif i == 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i, j + 1))
                    # append!(index_k, 1)

                    if epsilon12 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon12)
                    end

                    b[index_u(i, j)] = hx*left
                    b[index_v(i, j)] = 0

                elseif i == Nx
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i - 1, j))
                    append!(index_k, -1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i, j - 1))
                    # append!(index_k, 1)

                    if epsilon34 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon34)
                    end

                    b[index_u(i, j)] = hx*right
                    b[index_v(i, j)] = 0
                elseif j == 0
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -1)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, 1)

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, 1)

                    # append!(index_i, index_v(i, j))
                    # append!(index_j, index_v(i + 1, j))
                    # append!(index_k, 1)

                    if epsilon14 != 0
                        append!(index_i, index_u(i, j))
                        append!(index_j, index_u(i, j))
                        append!(index_k, epsilon14)
                    end

                    b[index_u(i, j)] = hx*below
                    b[index_v(i, j)] = 0
                else
                    # nonlinear terms
                    phi_x = x[index_v(i + 1, j)] - x[index_v(i - 1, j)]
                    phi_y = x[index_v(i, j + 1)] - x[index_v(i, j - 1)]

                    xi_x = x[index_u(i + 1, j)] - x[index_u(i - 1, j)]
                    xi_y = x[index_u(i, j + 1)] - x[index_u(i, j - 1)]

                    # right hand side
                    # first equation
                    b[index_u(i, j)] = x[index_u(i, j)] + tau / (2 * hx ^ 2) * (
                        x[index_u(i + 1, j)]
                        - 4 * x[index_u(i, j)]
                        + x[index_u(i - 1, j)]
                        + x[index_u(i, j + 1)]
                        + x[index_u(i, j - 1)]
                    )

                    # second equation
                    b[index_v(i, j)] = -(x[index_u(i, j)])
                    b[index_v(i, j)] -= (x[index_v(i + 1, j)]) / (hx ^ 2)
                    b[index_v(i, j)] -= (x[index_v(i - 1, j)]) / (hx ^ 2)
                    b[index_v(i, j)] -= (x[index_v(i, j + 1)]) / (hx ^ 2)
                    b[index_v(i, j)] -= (x[index_v(i, j - 1)]) / (hx ^ 2)
                    b[index_v(i, j)] += 4 * (x[index_v(i, j)]) / (hx ^ 2)

                    # diagonal of matrix
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, 1 + 2.0 * tau / (hx ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j))
                    append!(index_k, -4 / (hx ^ 2))

                    # u x
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j + 1))
                    append!(index_k, tau * Re / (8 * hx * hy) * phi_x - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, -tau * Re / (8 * hx * hy) * xi_x)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j - 1))
                    append!(index_k, -tau * Re / (8 * hx * hy) * phi_x - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j - 1))
                    append!(index_k, tau * Re / (8 * hx * hy) * xi_x)

                    # u y
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i + 1, j))
                    append!(index_k, -tau * Re / (8 * hx * hy) * phi_y - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, tau * Re / (8 * hx * hy) * xi_y)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i - 1, j))
                    append!(index_k, tau * Re / (8 * hx * hy) * phi_y - tau / (2.0 * hx ^ 2))

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i - 1, j))
                    append!(index_k, -tau * Re / (8 * hx * hy) * xi_y)

                    # equation v
                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, 1 / (hy ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i, j - 1))
                    append!(index_k, 1 / (hy ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, 1 / (hx ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_v(i - 1, j))
                    append!(index_k, 1 / (hx ^ 2))

                    append!(index_i, index_v(i, j))
                    append!(index_j, index_u(i, j))
                    append!(index_k, 1)
                end
            end
        end

        return (sparse(index_i, index_j, index_k, 2 * (Nx + 1) * (Ny + 1), 2 * (Nx + 1) * (Ny + 1)), b)
    end

    # start
    start_time = Dates.now()

    # y is an initial solution
    if isnothing(initial)
        y = zeros(2 * (Nx + 1) * (Ny + 1))
    else
        y = initial
        iteration = iteration + 1
    end

    # Loop updating x from x=A^-1*b
    # x is solution at time step n+1
    # y is solution at time step n
    x = deepcopy(y)
    error = 1
    # iteration = 1
    error_save = []
    while error > 1e-6 && iteration <= max_iteration
        if disp_error
            print("iteration: {}".format(iteration))
        end

        if iteration == 1
            A, b = create_matrix_stream(x, x, Nx=Nx, Ny=Ny, tau=tau, Re=Re, epsilon=epsilon)
            luA = lu(A)
            x = luA \ b

            @info "$(Dates.format(Dates.now(), Dates.dateformat"dd u Y -- H:M:S")) iter: $iteration/$max_iteration"
            println("\t max \t abs u: $(maximum(x[1:2:end]))")
            println("\t max \t abs v: $(maximum(x[2:2:end]))")
            iteration += 1
        else
            A, b = create_matrix_stream(x, y, Nx=Nx, Ny=Ny, tau=tau, Re=Re, epsilon=epsilon)

            # save time step n to y
            y = deepcopy(x)

            # solve linear sys
            luA = lu(A)
            x = luA \ b

            # calculate error
            @info "$(Dates.format(Dates.now(), Dates.dateformat"dd u Y -- H:M:S")) iter: $iteration/$max_iteration"
            println("\t max \t abs u: $(maximum(x[1:2:end]))")
            println("\t max \t abs v: $(maximum(x[2:2:end]))")

            u = x[1:2:end]
            v = x[2:2:end]
            error = maximum(abs.(x[2:2:end]-y[2:2:end]))
            println("\t error \t abs v: $(error)")
            push!(error_save, error)

            # save every 10 iterations
            if iteration % save_every == 0
                # plot_contour(reshape(v, (Nx + 1, Ny + 1)), "fig/save/checkpoint_$iteration.png")
                save_checkpoint(joinpath(location, "iter_u_$(iteration-1).txt"), reshape(u, (Nx + 1, Ny + 1)))
                save_checkpoint(joinpath(location, "iter_v_$(iteration-1).txt"), reshape(v, (Nx + 1, Ny + 1)))
            end
            
            iteration += 1
        end
    end
    
    u = x[1:2:end]
    v = x[2:2:end]
    u = reshape(u, (Nx + 1, Ny + 1))
    v = reshape(v, (Nx + 1, Ny + 1))
    save_checkpoint(joinpath(location, "iter_u_$iteration.txt"), u)
    save_checkpoint(joinpath(location, "iter_v_$iteration.txt"), v)
    

    stop_time = Dates.now()
    elapsed_time = Dates.canonicalize(Dates.CompoundPeriod(Dates.DateTime(stop_time) - Dates.DateTime(start_time)))
    @info "elapsed time = $elapsed_time"

    return u, v, elapsed_time, error_save, iteration
end


function plot_contour(v::Matrix; save_name=nothing, title=nothing)
    m = minimum(v)
    M = maximum(v)

    sc1 = collect(m:-m/20:0)
    sc2 = collect(0:M/100:M)
    sc = [sc1; sc2]

    # plot contour
    # Plots.savefig(Plots.contour(u, fill=false, levels=sc), "test_u.pdf");
    plt = Plots.contour(v, fill=false, levels=sc, color=:turbo), save_name

    # add title
    if !isnothing(title)
        title!(title)
    end

    # save 
    if !isnothing(save_name)
        Plots.savefig(plt, save_name)
    end

end


function save_checkpoint(file_name::String, matrix)
    open(file_name, "w") do io
        writedlm(io, matrix)
    end
end


function load_solution_txt(file_name::String)
    readdlm(file_name)
end


"""
    natural(x, y)

sort String on natural way e.g. A2 < A10

# Examples:
```julia-repl
julia> sort(["a1", "a2", "a10"], lt=VRPTW.natural)
3-element Vector{String}:
 "a1"
 "a2"
 "a10"
```
"""
function natural(x, y)
    k(x) = [occursin(r"\d+", s) ? parse(Int, s) : s
            for s in split(replace(x, r"\d+" => s -> " $s "))]
    A = k(x)
    B = k(y)
    for (a, b) in zip(A, B)
        if !isequal(a, b)
            return typeof(a) <: typeof(b) ? isless(a, b) :
                   isa(a, Int) ? true : false
        end
    end
    return length(A) < length(B)
end


function natural_sort(i)
    sort(i, lt=natural)
end


function plot_gif(; save_name="anim_fps5.gif", folder="checkpoints")
    # file_names = ["checkpoint/iter_$iter.txt" for iter in 5:5:100]
    file_names = sort(glob("*", folder), lt=natural)
    anim = @animate for file_name in file_names
        v = load_solution_txt(file_name)
        plot_contour(v, title=file_name)
    end
    gif(anim, save_name, fps = 5)
end