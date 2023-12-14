import REPL
using REPL.TerminalMenus



function start()

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

