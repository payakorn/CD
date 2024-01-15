import REPL
using REPL.TerminalMenus



function start1()

    options = [
        "Δx = Δy = τ = 0.1  and Re=100",
        "Δx = Δy = τ = 0.01 and Re=100",
        "Input Manually",
    ]
    
    
    menu = RadioMenu(options, pagesize=10)
    choice = request("Choose your favorite startup:", menu)

    print(" > Input max iteration (defalse=100): ")
    max_iteration = try readline() catch ArgumentError; nothing end
    if max_iteration |> isempty
        max_iteration = 100
    else
        max_iteration = parse(Int64, max_iteration)
    end
    
    if choice != -1
        if choice == 1
            hx=hy=tau=0.1
            Re=100

            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration);
        elseif choice == 2
            hx=hy=tau=0.01
            Re=100

            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration);
        elseif choice == 3
            print(" > Input Δx: ")
            hx = parse(Float64, readline())
            
            print(" > Input Δy: ")
            hy = parse(Float64, readline())

            print(" > Input  τ: ")
            tau = parse(Float64, readline())

            print(" > Input Re: ")
            Re = parse(Float64, readline())
            
            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration);
        end
    else
        println("Canceled.")
    end
    return u, v, elapsed_time, error_save, iteration
end


function start2()

    options = [
        "Δx = Δy = τ = 0.1  and Re=100",
        "Δx = Δy = τ = 0.01 and Re=100",
        "Δx = Δy = τ = 0.01 and Re=1000",
        "Input Manually",
    ]
    
    
    menu = RadioMenu(options, pagesize=10)
    choice = request("Choose your favorite startup:", menu)
    
    
    # 
    options2 = [
        "Each epsilon = 1e-8",
        "Input the same epsilon:",
        "Input Manually",
    ]
    menu2 = RadioMenu(options2, pagesize=10)
    choice2 = request("Choose your favorite epsilon:", menu2)

    if choice2 == 1
        epsilon1  = 1e-8
        epsilon2  = 1e-8
        epsilon3  = 1e-8
        epsilon4  = 1e-8
        epsilon12 = 1e-8
        epsilon23 = 1e-8
        epsilon34 = 1e-8
        epsilon14 = 1e-8

    elseif choice2 == 2
        print(" > Input all epsilon: ")
        epsilon1 = parse(Float64, readline())
        epsilon2  = epsilon1
        epsilon3  = epsilon1
        epsilon4  = epsilon1
        epsilon12 = epsilon1
        epsilon23 = epsilon1
        epsilon34 = epsilon1
        epsilon14 = epsilon1

    elseif choice2 == 3
        print(" > Input epsilon1: ")
        epsilon1 = parse(Float64, readline())
        
        print(" > Input epsilon2: ")
        epsilon2 = parse(Float64, readline())

        print(" > Input epsilon3: ")
        epsilon3 = parse(Float64, readline())

        print(" > Input epsilon4: ")
        epsilon4 = parse(Float64, readline())

        print(" > Input epsilon12: ")
        epsilon12 = parse(Float64, readline())

        print(" > Input epsilon23: ")
        epsilon23 = parse(Float64, readline())

        print(" > Input epsilon34: ")
        epsilon34 = parse(Float64, readline())

        print(" > Input epsilon14: ")
        epsilon14 = parse(Float64, readline())

    end

    epsilon = (
        epsilon1,
        epsilon2,
        epsilon3,
        epsilon4,
        epsilon12,
        epsilon23,
        epsilon34,
        epsilon14,
    )


    print(" > Input max iteration (default=100): ")
    max_iteration = try readline() catch ArgumentError; nothing end
    if max_iteration |> isempty
        max_iteration = 100
    else
        max_iteration = parse(Int64, max_iteration)
    end
    
    if choice != -1
        if choice == 1
            hx=hy=tau=0.1
            Re=100

            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration, epsilon=epsilon);
        elseif choice == 2
            hx=hy=tau=0.01
            Re=100

            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration, epsilon=epsilon);
        elseif choice == 3
            hx=hy=tau=0.01
            Re=1000

            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration, epsilon=epsilon);
        elseif choice == 4
            print(" > Input Δx: ")
            hx = parse(Float64, readline())
            
            print(" > Input Δy: ")
            hy = parse(Float64, readline())

            print(" > Input  τ: ")
            tau = parse(Float64, readline())

            print(" > Input Re: ")
            Re = parse(Float64, readline())
            
            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration, epsilon=epsilon);
        end
    else
        println("Canceled.")
    end
    return u, v, elapsed_time, error_save, iteration-1
end


function seleting_start()

    """
    
        first selecting model
    
    """

    options_start = [
        "6_epsilon",
        "8_epsilon",
    ]
    menu_start = RadioMenu(options_start, pagesize=10)
    choice_start = request("Choose your favorite start:", menu_start)


    if choice_start == 1
        @info "6 epsilon"
        sys_stream = sys_stream1

        options2 = [
        "Each epsilon = 1e-8",
        "Input the same epsilon:",
        "Input Manually",
        ]
        menu2 = RadioMenu(options2, pagesize=10)
        choice2 = request("Choose your favorite epsilon:", menu2)

        if choice2 == 1

            epsilon1  = 1e-8
            epsilon2  = 1e-8
            epsilon3  = 1e-8
            epsilon4  = 1e-8
            epsilon5 = 1e-8
            epsilon6 = 1e-8
            
        elseif choice2 == 2

            print(" > Input all epsilon: ")
            epsilon1 = parse(Float64, readline())
            epsilon2  = epsilon1
            epsilon3  = epsilon1
            epsilon4  = epsilon1
            epsilon5 = epsilon1
            epsilon6 = epsilon1

        elseif choice2 == 3

            print(" > Input epsilon1: ")
            epsilon1 = parse(Float64, readline())
            
            print(" > Input epsilon2: ")
            epsilon2 = parse(Float64, readline())

            print(" > Input epsilon3: ")
            epsilon3 = parse(Float64, readline())

            print(" > Input epsilon4: ")
            epsilon4 = parse(Float64, readline())

            print(" > Input epsilon5: ")
            epsilon5 = parse(Float64, readline())

            print(" > Input epsilon6: ")
            epsilon6 = parse(Float64, readline())

        end

        epsilon = (
            epsilon1,
            epsilon2,
            epsilon3,
            epsilon4,
            epsilon5,
            epsilon6,
        )

    elseif choice_start == 2

        @info "8 epsilon"
        sys_stream = sys_stream2

        options2 = [
            "Each epsilon = 1e-8",
            "Input the same epsilon:",
            "Input Manually",
        ]
        menu2 = RadioMenu(options2, pagesize=10)
        choice2 = request("Choose your favorite epsilon:", menu2)

        if choice2 == 1

            epsilon1  = 1e-8
            epsilon2  = 1e-8
            epsilon3  = 1e-8
            epsilon4  = 1e-8
            epsilon12 = 1e-8
            epsilon23 = 1e-8
            epsilon34 = 1e-8
            epsilon14 = 1e-8

        elseif choice2 == 2

            print(" > Input all epsilon: ")
            epsilon1 = parse(Float64, readline())
            epsilon2  = epsilon1
            epsilon3  = epsilon1
            epsilon4  = epsilon1
            epsilon12 = epsilon1
            epsilon23 = epsilon1
            epsilon34 = epsilon1
            epsilon14 = epsilon1

        elseif choice2 == 3

            print(" > Input epsilon1: ")
            epsilon1 = parse(Float64, readline())
            
            print(" > Input epsilon2: ")
            epsilon2 = parse(Float64, readline())

            print(" > Input epsilon3: ")
            epsilon3 = parse(Float64, readline())

            print(" > Input epsilon4: ")
            epsilon4 = parse(Float64, readline())

            print(" > Input epsilon12: ")
            epsilon12 = parse(Float64, readline())

            print(" > Input epsilon23: ")
            epsilon23 = parse(Float64, readline())

            print(" > Input epsilon34: ")
            epsilon34 = parse(Float64, readline())

            print(" > Input epsilon14: ")
            epsilon14 = parse(Float64, readline())

        end

        epsilon = (
            epsilon1,
            epsilon2,
            epsilon3,
            epsilon4,
            epsilon12,
            epsilon23,
            epsilon34,
            epsilon14,
        )
    end

    return options_start[choice_start], sys_stream, epsilon

end


function selecting_settings(;epsilon=nothing, mode=nothing)

    """
    
        second section (selecting => dx, dy, tau, Re, epsilon)
    
    """

    options = [
        "Δx = Δy = 0.1 , τ = Δx² and Re=100",
        "Δx = Δy = 0.01, τ = Δx² and Re=100",
        "Δx = Δy = 0.01, τ = Δx² and Re=1000",
        "Input Manually",
    ]
    menu = RadioMenu(options, pagesize=10)
    choice = request("Choose your favorite startup:", menu)
    
    # select setting choice
    if choice == 1

        hx = hy = 0.1
        tau = hx^2
        Re = 100

    elseif choice == 2

        hx = hy = 0.01
        tau = hx^2
        Re = 100

    elseif choice == 3

        hx = hy = 0.01
        tau = hx^2
        Re = 1000

    elseif choice == 4

        print(" > Input Δx: ")
        hx = parse(Float64, readline())
        
        print(" > Input Δy: ")
        hy = parse(Float64, readline())

        print(" > Input  τ: ")
        tau = parse(Float64, readline())

        print(" > Input Re: ")
        Re = parse(Float64, readline())
        
    end

    if isnothing(epsilon) || isnothing(mode)
        print(" > Input the maximum number of iterations (default=100): ")
    else
        print(" > Input the maximum number of iterations (default=100) | currently have $(count_save_history(mode, hx, hy, tau, Re, epsilon...)) iterations: ")
    end
    max_iteration = try readline() catch ArgumentError; nothing end

    if max_iteration |> isempty
        max_iteration = 100
    else
        max_iteration = parse(Int64, max_iteration)
    end

    return hx, hy, tau, Re, max_iteration

end


function get_latest_u(mode::String, Δx::Float64, Δy::Float64, τ::Float64, Re::Real, epsilon...)
    location = folder(mode, Δx, Δy, τ, Re, epsilon...)
    all_files = glob("iter_u_*.txt", location) |> natural_sort
    return load_solution_txt(all_files[end])
end


function get_latest_v(mode::String, Δx::Float64, Δy::Float64, τ::Float64, Re::Real, epsilon...)
    location = folder(mode, Δx, Δy, τ, Re, epsilon...)
    all_files = glob("iter_v_*.txt", location) |> natural_sort
    return load_solution_txt(all_files[end])
end


function get_latest_uv(mode::String, Δx::Float64, Δy::Float64, τ::Float64, Re::Real, epsilon...)

    # load matrix u and v --> convert to vector
    u = get_latest_u(mode, Δx, Δy, τ, Re, epsilon...) |> vec
    v = get_latest_v(mode, Δx, Δy, τ, Re, epsilon...) |> vec

    x = zeros(length(u) + length(v))

    x[1:2:end] = u
    x[2:2:end] = v

    return x

end


function start()
    
    
    starting_choice, sys_stream, epsilon = seleting_start()
    hx, hy, tau, Re, max_iteration = selecting_settings(epsilon=epsilon, mode=starting_choice)

    # check if the max_iteration > the previous run
    latest_iteration = count_save_history(starting_choice, hx, hy, tau, Re, epsilon...)
    if latest_iteration >= max_iteration
        @error "this setting is already completed (currently have >= the maximum number of iterations) => try again"
        return nothing
    else

        if latest_iteration == 0

            x = nothing
            iteration = 1

        else

            options = [
                "run using the latest save solution as initial solution",
                "run from zero initial solution",
            ]
            menu = RadioMenu(options, pagesize = 10)
            initial_choice = request("Choose your favorite options: ", menu)

            if initial_choice == 1

                x = get_latest_uv(starting_choice, hx, hy, tau, Re, epsilon...)
                iteration = latest_iteration

            elseif initial_choice == 2
                
                x = nothing
                iteration = 1

            end

        end

    end

    # solving problem
    u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration, epsilon=epsilon, initial=x, iteration=iteration)

    # save history
    save_history(starting_choice, hx, hy, tau, Re, elapsed_time, iteration-1)

    # plot  error
    error_plot = true
    if error_plot && length(error_save) > 1
        
        len = length(error_save)
        in_x = fld(len, 2):10:len
        plt = scatter(in_x, error_save[in_x])
        title!("$(starting_choice)_x=$(hx)_y=$(hy)_tau=$(tau)_Re=$(Re)_eps=$(epsilon[1])")
        Plots.savefig(plt, "fig/$(starting_choice)_x=$(hx)_y=$(hy)_tau=$(tau)_Re=$(Re)_eps=$(epsilon[1])_iter=$(iteration-1).png")

    end

    return u, v, elapsed_time, error_save, iteration-1
end


function save_history(mode::String, Δx::Float64, Δy::Float64, τ::Float64, Re::Real, elapsed_time::Dates.CompoundPeriod, iteration::Integer)

    # change format of elapsed_time
    sti = replace("$(elapsed_time)", ' ' => '-')
    elapsed_time = replace(sti, ','=> '-')

    head_name = "date,mode,hx,hy,tau,Re,elapsed_time,iteration\n"
    if !isfile("history.csv")

        open("history.csv", "w") do file
            write(file, head_name)
            data = "$(now()),$mode,$Δx,$Δy,$τ,$Re,$elapsed_time,$iteration\n"
            write(file, data)
        end

    else

        open("history.csv", "a") do file
            data = "$(now()),$mode,$Δx,$Δy,$τ,$Re,$elapsed_time,$iteration\n"
            write(file, data)
        end

    end

end


function folder(mode::String, Δx::Float64, Δy::Float64, τ::Float64, Re::Real, epsilon...)
    return joinpath("checkpoints", mode, "$(Δx)_$(Δy)_$(τ)_$(Re)_$(epsilon[1])")
end



"""
    count_save_history(mode::String, Δx::Float64, Δy::Float64, τ::Float64, Re::Real, epsilon...)

    return the last number of iteration in the folder
"""
function count_save_history(mode::String, Δx::Float64, Δy::Float64, τ::Float64, Re::Real, epsilon...)

    """
        1. find all files in the folder and sort them to get the latest file
        2. split directory to obtain the files name (file name in the form iter_u_[iter].txt)
        3. get the number of [iter]
    """

    location = folder(mode, Δx, Δy, τ, Re, epsilon...)
    all_files = glob("iter_u_*", location) |> natural_sort
    if isempty(all_files)
        return 0
    else
        the_last_file = splitdir(all_files[end])[end]
        the_number_of_last_file = parse(Int64, split(split(the_last_file, ".")[1], "_")[end])

        return the_number_of_last_file
    end
end
