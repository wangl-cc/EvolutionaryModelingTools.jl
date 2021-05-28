using GillespieTools
using RecordedArrays

@args birth_c(r, x) = r * x
@args birth_u!(ind, x) = x[ind] += 1

@args comp_c(r, K, x) = r * x * x' / K
@args comp_u!(ind, x) = x[ind[1]] -=1

birth = Reaction(birth_c, birth_u!)
comp = Reaction(comp_c, comp_u!)

c = ContinuousClock(100.0)
ps = (r = 0.5, K = 100, x = DynamicRArray(c, 10))
rs = (birth, comp)

gillespie(c, ps, rs)
