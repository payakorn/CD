import REPL
using REPL.TerminalMenus



function start()

    options = [
        "hx=hy=τ=0.1, Re=100, maxiter=10",
        "hx=hy=τ=0.01, Re=100, maxiter=100",
        "Manually Input"
    ]
    
    
    menu = RadioMenu(options, pagesize=10)
    choice = request("Choose your favorite startup:", menu)
    
    if choice != -1
        if choice == 1
            hx=hy=tau=0.1
            Re=100
            maxiter=10
            solve(hx, hy, tau, Re; maxiter=maxiter);
        elseif choice == 2
            hx=hy=tau=0.01
            Re=100
            maxiter=100
            solve(hx, hy, tau, Re; maxiter=maxiter);
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
            maxiter = parse(Int64, readline())

            solve(hx, hy, tau, Re; maxiter=maxiter);
        end
    else
        println("Canceled.")
    end

    # hx, hy, tau, Re; maxiter=10
end

