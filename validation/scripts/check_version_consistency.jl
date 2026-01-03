using OpenPKPDCore

root = normpath(joinpath(@__DIR__, ".."))
version_file = joinpath(root, "VERSION")

v = strip(read(version_file, String))

if v != OpenPKPDCore.OPENPKPD_VERSION
    error("VERSION file ($(v)) does not match OpenPKPDCore.OPENPKPD_VERSION ($(OpenPKPDCore.OPENPKPD_VERSION))")
end

println("Version consistency check passed: " * v)
