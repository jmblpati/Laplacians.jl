### samplingSolver
include("revampedLinkedListFloatStorage.jl")
export llsInit
export llsAdd
export llsPurge

include("sqLinOpWrapper.jl")
export SqLinOp

include("samplingSolver.jl")
#export samplingSDDSolver
#export samplingLapSolver
#export samplingParams

### KMP
include("KMPSolver.jl")
export KMPSDDSolver
export KMPLapSolver
export KMPParams


### useful misc
include("condNumber.jl")
export condNumber

include("powerIteration.jl")
export powerIteration

"""A list containing SDD linear system solvers. They take in a SDD matrix plus tol, maxits and maxtime parameters."""
SDDSolvers = [augTreeSolver, KMPSDDSolver, AMGSolver]

"""A list containing Laplacian linear system solvers. They take in an adjacency matrix plus tol, maxits and maxtime parameters."""
LapSolvers = [augTreeLapSolver, KMPLapSolver, AMGLapSolver]
