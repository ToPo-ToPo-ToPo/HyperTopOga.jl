
#
include("base_kit/mmasub.jl")
include("base_kit/subsolv.jl")
include("base_kit/kktcheck.jl")

# GCMMA only
include("base_kit/gcmmasub.jl")
include("base_kit/asymp.jl")
include("base_kit/raaupdate.jl")
include("base_kit/concheck.jl")

#
include("history.jl")

# 共有部分
include("abstract_MMA/abstract_MMA.jl")
include("abstract_MMA/initial_compute.jl")
include("abstract_MMA/compute.jl")
include("abstract_MMA/output_vtu_element_base.jl")
include("abstract_MMA/output_vtu.jl")

# GCMMA
include("GCMMA_wrapper.jl")