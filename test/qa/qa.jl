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
            # Names accessed qualified that are still un-`public`/un-exported upstream.
            # (SciMLBase-owned names are now public in 3.27.0 and need no ignore;
            #  `Random.GLOBAL_RNG` was migrated to the public `Random.default_rng()`.)
            # Each remaining name is a genuine external API that upstream has not yet
            # declared `public`; verified non-public against the resolved releases
            # (Symbolics 7.29.0, SciMLStructures 1.10.1) on Julia 1.12.
            ignore = (
                :StanTarget,             # Symbolics build-target type, only `import`ed
                                         # (not exported/`public`) into ModelingToolkit
                :Tunable,                # SciMLStructures core API, not yet `public`
                :canonicalize,           # SciMLStructures core API, not yet `public`
                :isscimlstructure,       # SciMLStructures core API, not yet `public`
            ),
        ),
    ),
    # DiffEqBayes pulls 68 names implicitly from heavy `using` deps (Turing,
    # Distributions, ModelingToolkit, DiffEqBase, ...); making them explicit is a
    # large, risky refactor tracked in https://github.com/SciML/DiffEqBayes.jl/issues/391
    ei_broken = (:no_implicit_imports,)
)
