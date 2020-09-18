"""
package for Variational mode decomposition
"""
module VMD
    using FFTW
    using Printf
    using RecipesBase
    
    include("svmd.jl")
    include("mvmd.jl")

end
