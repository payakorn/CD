import REPL
using REPL.TerminalMenus



function start()

    options = [
        "hx=hy=τ=0.1, Re=100, maxiter=100",
        "hx=hy=τ=0.01, Re=100, maxiter=100",
        "Manually Input"
    ]
    
    
    menu = RadioMenu(options, pagesize=10)
    choice = request("Choose your favorite startup:", menu)
    
    if choice != -1
        if choice == 1
            hx=hy=tau=0.1
            Re=100
            max_iteration=100
            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration);
        elseif choice == 2
            hx=hy=tau=0.01
            Re=100
            max_iteration=100
            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration);
        elseif choice == 3
            print("Input hx: ")
            hx = parse(Float64, readline())
            
            print("Input hy: ")
            hy = parse(Float64, readline())

            print("Input τ: ")
            tau = parse(Float64, readline())

            print("Input Re: ")
            Re = parse(Float64, readline())
            
            print("Input max iter: ")
            max_iteration = parse(Int64, readline())

            u, v, elapsed_time, error_save, iteration = sys_stream(hx, hy, tau, Re; max_iteration=max_iteration);
        end
    else
        println("Canceled.")
    end
    return u, v, elapsed_time, error_save, iteration
end

