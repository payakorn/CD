# new section for Stream function
function sys_stream(hx, hy, tau, Re; max_iteration=100, disp_error=false, epsilon=false, left=0, right=0, top=1, below=0)
    # input step size (hx, hy) and (tau)
    # save = True if we need to export solution

    # start_time = time.time()

    xr = 0
    # xr -= hx/2
    xl = 1
    # xl += hx/2

    yr = 0
    # yr -= hy/2
    yl = 1
    # yl += hy/2

    # the number
    Nx = Int64((xl - xr) / (hx))
    Ny = Int64((yl - yr) / (hy))
    # if disp_error is True:
    #     print("Nx = {}".format(Nx))
    #     print("Ny = {}".format(Ny))

    # function return the position of u and in vector
    function index_u(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 1
    end

    function index_v(i, j)
        return 2 * ((i) * (Ny + 1) + j) + 2
    end

    # create empty vector for storing element in matrix
    index_i = []
    index_j = []

    function create_matrix_stream(x, y; Nx=false, Ny=false, tau=false, Re=Re, index_i=index_i, index_j=index_j)
        index_i = []
        index_j = []
        index_k = []

        # create empty right hamd side vector
        b = zeros(2 * (Nx + 1) * (Ny + 1))

        # epsilon
        if !epsilon
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

                    # add
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j + 1))
                    append!(index_k, -tau * Re / (8 * hx * hy) * xi_x)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i, j - 1))
                    append!(index_k, -tau * Re / (8 * hx * hy) * phi_x - tau / (2.0 * hx ^ 2))

                    # add
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i, j - 1))
                    append!(index_k, tau * Re / (8 * hx * hy) * xi_x)

                    # u y
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i + 1, j))
                    append!(index_k, -tau * Re / (8 * hx * hy) * phi_y - tau / (2.0 * hx ^ 2))

                    # add
                    append!(index_i, index_u(i, j))
                    append!(index_j, index_v(i + 1, j))
                    append!(index_k, tau * Re / (8 * hx * hy) * xi_y)

                    append!(index_i, index_u(i, j))
                    append!(index_j, index_u(i - 1, j))
                    append!(index_k, tau * Re / (8 * hx * hy) * phi_y - tau / (2.0 * hx ^ 2))

                    # add
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

    # y is initial solution
    y = zeros(2 * (Nx + 1) * (Ny + 1))

    # Loop updating x from x=A^-1*b
    # x is solution at time step n+1
    # y is solution at time step n
    x = deepcopy(y)
    error = 1
    iteration = 1
    error_save = []
    while error > 1e-6 && iteration <= max_iteration
        if disp_error
            print("iteration: {}".format(iteration))
        end

        if iteration == 1
            A, b = create_matrix_stream(x, x, Nx=Nx, Ny=Ny, tau=tau, Re=Re)
            luA = lu(A)
            x = luA \ b
            # error = np.linalg.norm(x - y)
            iteration += 1
        else
            A, b = create_matrix_stream(x, y, Nx=Nx, Ny=Ny, tau=tau, Re=Re)

            # save time step n to y
            y = deepcopy(x)

            # solve linear sys
            luA = lu(A)
            x = luA \ b

            # calculate error
            # u = [x[i] for i in range(0, 2 * (Nx + 1) * (Ny + 1), 2)]
            # v = [x[i] for i in range(1, 2 * (Nx + 1) * (Ny + 1), 2)]
            # error = np.amax(np.abs([x[i] - y[i] for i in range(1, 2 * (Nx + 1) * (Ny + 1), 2)]))
            # error_save.append(error)

            # if disp_error
            #     print("max   abs u: {}".format(np.amax(np.abs(u))))
            #     print("max   abs v: {}".format(np.amax(np.abs(v))))
            #     print("error abs v: {}".format(error))
            # end
            @info "iter: $iteration/$max_iteration"
            iteration += 1
        end
    end

    u = x[2:2:end]
    v = x[1:2:end]

    u = reshape(u, (Nx + 1, Ny + 1))
    v = reshape(v, (Nx + 1, Ny + 1))

    # elapsed time
    # elapsed_time = time.time() - start_time
    # print("time = ", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    # u is even elements of x
    # u = [x[i] for i in range(0, 2 * (Nx + 1) * (Ny + 1), 2)]
    # u = np.reshape(u, [Nx + 1, Ny + 1])

    # # v is odd elements of x
    # v = [x[i] for i in range(1, 2 * (Nx + 1) * (Ny + 1), 2)]
    # v = np.reshape(v, [Nx + 1, Ny + 1])

    # # this for export to CSV
    # if save
    #     (X, Y) = grid(hx, hy, r=1)
    #     np.savetxt("X_stream_{}_{}_{}.csv".format(Nx, Ny, Re), X, delimiteration=",")
    #     np.savetxt("Y_stream_{}_{}_{}.csv".format(Nx, Ny, Re), Y, delimiteration=",")
    #     np.savetxt("xi_{}_{}_{}.csv".format(Nx, Ny, Re), u, delimiteration=",")
    #     np.savetxt("phi_{}_{}_{}.csv".format(Nx, Ny, Re), v, delimiteration=",")
    # end

    elapsed_time = 0

    # return u, v, elapsed_time, error_save, iteration
    return u, v, elapsed_time, error_save, iteration
end


# def stream( hx, hy, tau, Re; save=false, max_iteration=100, disp_error=false, epsilon=false, left=0, right=0, top=1, below=0)
    # u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re, max_iteration=max_iteration, disp_error=disp_error, epsilon=epsilon, left=left, right=right, top=top, below=below)
    # max_u = np.amax(np.abs(u))
    # max_v = np.amax(np.abs(v))

    # # save output to text file
    # if save is True:
    #     np.savetxt(
    #         "stream_save/u-{}-{}-{}.txt".format(Re, hx, iteration), u, fmt="%1.4e"
    #     )  # use exponential notation
    #     np.savetxt(
    #         "stream_save/v-{}-{}-{}.txt".format(Re, hx, iteration), v, fmt="%1.4e"
    #     )  # use exponential notation
    #     np.savetxt(
    #         "stream_save/u-{}-{}-{}.out".format(Re, hx, iteration), u, fmt="%1.4e"
    #     )  # use exponential notation
    #     np.savetxt(
    #         "stream_save/v-{}-{}-{}.out".format(Re, hx, iteration), v, fmt="%1.4e"
    #     )  # use exponential notation
    #     np.savetxt(
    #         "stream_save/v-error-{}-{}.txt".format(Re, hx), error_save, fmt="%1.4e"
    #     )  # use exponential notation
    #     np.savetxt(
    #         "stream_save/v-error-{}-{}.out".format(Re, hx), error_save, fmt="%1.4e"
    #     )  # use exponential notation

    # if plot is True:
    #     # X, Y = grid(hx, hy, r=1)
    #     xi = np.transpose(u)
    #     phi = np.transpose(v)
    #     m = np.amin(phi)
    #     M = np.amax(phi)
    #     print("min phi = {}".format(m))
    #     print("max phi = {}".format(M))
    #     # cp = plt.contourf(x, y, phi, cmap='RdGy', )
    #     sort_phi = np.sort(np.reshape(phi, -1))
    #     N = np.argmax(sort_phi > 0)
    #     new_sc1 = sort_phi[0:N-500:100]
    #     new_sc2 = sort_phi[N:-1]
    #     sc1 = np.arange(m, 0, -m / 20)
    #     sc2 = np.arange(0, M, M / 20)
    #     new_sc = np.append(new_sc1, new_sc2)
    #     sc = np.append(sc1, sc2)
    #     cp = plt.contour(phi, levels=sc, cmap="RdGy")
    #     plt.gca().set_aspect("equal")
    #     print("fig size = {}".format(np.shape(phi)))
    #     plt.clabel(cp, inline=True, fontsize=8)
    #     plt.colorbar(cp)
    #     plt.show()
    #     # plot_3d(x, y, phi)

    # print("max norm xi  = {}".format(max_u))
    # print("max norm phi = {}".format(max_v))
    # print("error = {}".format(error_save[-1]))

    # return xi, phi, error_save

# def four_points(hx):
#     epsilon = (1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8)
#     u, v, error = stream(
#         hx=hx,
#         hy=hx,
#         tau=hx,
#         Re=100,
#         plot=True,
#         save=False,
#         max_iteration=500,
#         disp_error=False,
#         epsilon=epsilon,
#     )

#     xr = 0
#     xr -= hx/2
#     xl = 1
#     xl += hx/2

#     yr = 0
#     yr -= hx/2
#     yl = 1
#     yl += hx/2

#     Nx = int((xl - xr) / (hx))
#     Ny = int((yl - yr) / (hx))

#     xi = np.zeros()

#     for i in range(1, Nx):
#         for j in range(1, Ny):
#             xi[i, j] = ()
    


# # # run code here
# if __name__ == "__main__":
#     hx = 1 / 52
#     tau = hx^2
#     epsilon = (1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 0)
#     u, v, error = stream(
#         hx=hx,
#         hy=hx,
#         tau=tau,
#         Re=100,
#         plot=True,
#         save=False,
#         max_iteration=100,
#         disp_error=True,
#         epsilon=epsilon,
#         left=0,
#         right=0,
#         top=4,
#         below=0,
#     )
#     # plot boundray

    