using SciMLTesting, DiffEqBayes, Test
using JET

run_qa(
    DiffEqBayes;
    explicit_imports = true,
    ei_kwargs = (;
        all_qualified_accesses_via_owners = (;
            # Re-exports surfaced through ModelingToolkit (Base.which points at the
            # defining package, not the module DiffEqBayes accesses them from):
            ignore = (
                :StanTarget,   # owner Symbolics, re-exported by ModelingToolkit
                :get_iv,       # owner ModelingToolkitBase, re-exported by ModelingToolkit
            ),
        ),
        all_qualified_accesses_are_public = (;
            # Non-public (un-`public`/un-exported) names of other packages that
            # DiffEqBayes accesses qualified; safe until each upstream marks them public:
            ignore = (
                :AbstractSciMLProblem,   # SciMLBase
                :Fix1,                   # Base
                :GLOBAL_RNG,             # Random
                :StanTarget,             # ModelingToolkit (owner Symbolics)
                :Tunable,                # SciMLStructures
                :canonicalize,           # SciMLStructures
                :get_iv,                 # ModelingToolkit (owner ModelingToolkitBase)
                :isscimlstructure,       # SciMLStructures
            ),
        ),
    ),
    # DiffEqBayes pulls 68 names implicitly from heavy `using` deps (Turing,
    # Distributions, ModelingToolkit, DiffEqBase, ...); making them explicit is a
    # large, risky refactor tracked in https://github.com/SciML/DiffEqBayes.jl/issues/391
    ei_broken = (:no_implicit_imports,)
)
