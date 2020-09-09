# FourBarLinkages

```julia
using FourBarLinkages

# coupler points in isotropic coordinates
coupler_points = [complex(0.8961867,-0.09802917),
        complex(1.2156535, -1.18749100),
        complex(1.5151435, -0.85449808),
        complex(1.6754775,  -0.48768058),
        complex(1.7138690,-0.30099232),
        complex(1.7215236,0.03269953),
        complex(1.6642029, 0.33241088),
        complex(1.4984171, 0.74435576),
        complex(1.3011834,  0.92153806)]

# You can compute a generic solution and store it in a file with the following line:
#       compute_generic_solutions(; filename=data/four_bar_start_solutions.jld2)

# compute four bar linkages for the couple points from the stored results
fourbars = four_bars(coupler_point)

# pick a fourbar
F = fourbars[3]
# let's animate this with Makie
animate(F, coupler_points)
# create endless loop (interrupt to stop)
animate(F, coupler_points; loop=true)
# save animation and hide axis
animate(F, coupler_points; filename="four-bar.gif", show_axis=true)
```


DynamicPolynomials = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HomotopyContinuation = "f213a82b-91d6-5c5d-acf7-10f1c761b327"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
MultivariatePolynomials = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
StaticPolynomials = "62e018b1-6e46-5407-a5a7-97d4fbcae734"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
