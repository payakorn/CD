try using CD catch e; import Pkg; Pkg.activate("."); Pkg.instantiate(); using CD end

u, v, elapsed_time, error_save, iteration = CD.start();

# using Plots

# # u, v, time, er, iteration = CD.sys_stream(0.01, 0.01, 0.01, 100, max_iteration=100)
# m = minimum(v)
# M = maximum(v)

# sc1 = collect(m:-m/20:0)
# sc2 = collect(0:m/50:M)
# sc = [sc1; sc2]
# Plots.savefig(Plots.contour(u, fill=false, levels=sc), "test_u.pdf");
# Plots.savefig(Plots.contour(v, fill=false, levels=sc, color=:turbo), "test_v.pdf");
