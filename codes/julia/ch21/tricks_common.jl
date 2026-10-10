# tricks_common.jl
#
# Shared utilities for Chapter 21 (Special Tricks for the Spectral Researcher).
#
# Provides:
#   * NAVY, CORAL, TEAL, PURPLE, GOLD, OLIVE colour constants;
#   * `apply_theme!()` to set CairoMakie theme to CMU Serif at the
#     standard textbook size;
#   * `output_dir(script_dir)` to resolve textbook/figures/ch21/julia/;
#   * `write_json(path, x)` for the etude --dump option.
#
# Author: Dr. Denys Dutykh (Khalifa University, Abu Dhabi, UAE)
# Part of "Computational Etudes: A Spectral Approach"

using CairoMakie
using Colors

const NAVY   = colorant"#142D6E"
const SKY    = colorant"#7896D2"
const CORAL  = colorant"#E74C3C"
const TEAL   = colorant"#16A085"
const PURPLE = colorant"#8E44AD"
const ORANGE = colorant"#E67E22"
const GOLD   = colorant"#D4A017"
const OLIVE  = colorant"#6B8E23"

function apply_theme!()
    set_theme!(Theme(fontsize = 11,
        fonts = (regular = "CMU Serif", bold = "CMU Serif Bold",
                 italic = "CMU Serif Italic"),
        Axis = (xlabelsize = 11, ylabelsize = 11, titlesize = 11,
                xticklabelsize = 10, yticklabelsize = 10,
                spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8),
        Legend = (labelsize = 9,)))
end

function output_dir(script_dir::AbstractString)
    p = script_dir
    while !isempty(p) && !isdir(joinpath(p, "textbook"))
        new_p = dirname(p)
        new_p == p && break
        p = new_p
    end
    out = joinpath(p, "textbook", "figures", "ch21", "julia")
    mkpath(out)
    return out
end

# ---- Minimal JSON writer (no external dependency) ----
_json_val(x::Bool) = x ? "true" : "false"
_json_val(x::Integer) = string(x)
_json_val(x::Real) = isfinite(x) ? string(Float64(x)) :
    (isnan(x) ? "NaN" : (x > 0 ? "Infinity" : "-Infinity"))
_json_val(x::Complex) = "{\"re\": " * _json_val(real(x)) *
                        ", \"im\": " * _json_val(imag(x)) * "}"
_json_val(x::AbstractString) = "\"" *
    replace(String(x), "\\" => "\\\\", "\"" => "\\\"") * "\""
_json_val(v::AbstractVector) = "[" * join(_json_val.(v), ", ") * "]"
_json_val(t::Tuple) = "[" * join(_json_val.(collect(t)), ", ") * "]"
_json_val(d::AbstractDict) = "{" *
    join([_json_val(string(k)) * ": " * _json_val(v) for (k, v) in d], ", ") *
    "}"

write_json(path, x) = open(io -> print(io, _json_val(x)), path, "w")
