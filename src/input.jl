import REPL
using REPL.TerminalMenus



function start()

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
        "All epsilon = 1e-8",
        "Input Manually",
    ]
    menu2 = RadioMenu(options2, pagesize=10)
    choice2 = request("Choose your favorite epsilon:", menu2)

    if choice2 == 1
        epsilon1 = 1e-8
        epsilon2 = 1e-8
        epsilon3 = 1e-8
        epsilon4 = 1e-8
        epsilon12 = 1e-8
        epsilon23 = 1e-8
        epsilon34 = 1e-8
        epsilon14 = 1e-8
    elseif choice == 2
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

