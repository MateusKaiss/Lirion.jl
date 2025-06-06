abstract type OPFModel end

struct DC <: OPFModel end
struct AC <: OPFModel end
struct SOCP <: OPFModel end

abstract type SubproblemAction end

struct SubproblemBuild <: SubproblemAction end
struct SubproblemUpdate <: SubproblemAction end
struct SubproblemGetResult <: SubproblemAction end
struct SubproblemNameElementsDimensions <: SubproblemAction end
struct SubproblemWriteResult <: SubproblemAction end
