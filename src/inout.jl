function input()

    # hx, hy, tau, Re; maxiter=10

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

    solve(hx, hy, tau, Re; maxiter=maxiter)
end

