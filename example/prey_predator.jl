using GillespieTools
using RecordedArrays

@args birth_c(r, x) = r * x
@args birth_u!(ind, x) = x[ind] += 1

@args comp_c(r, K, x) = r * x * x' / K
@args comp_u!(ind, x) = x[ind[1]] -=1

@args pred_c(p, x, y) = p * x * y'
@args function pred_u!(ind, x, y)
    i, j = ind.I
    x[i] -= 1
    y[j] += 1
    nothing
end

@args death_c(d, y) = d * y
@args death_u!(ind, y) = y[ind] -= 1

c = ContinuousClock(100.0)

birth = Reaction(birth_c, birth_u!)
comp  = Reaction(comp_c, comp_u!)
pred  = Reaction(pred_c, pred_u!)
death = Reaction(death_c, death_u!)

x, y = DynamicRArray(c, [100], [10])
ps = (r=1.0, K=100, p=0.001, d=0.1, x, y)
rs = (birth, comp, pred, death)

gillespie(c, ps, rs)
