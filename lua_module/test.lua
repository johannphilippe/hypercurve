package.cpath = package.cpath .. ";/home/johann/Documents/GitHyb/build-hypercurve-Clang10-Debug"

local hc =  require("liblua_hypercurve")

local div = 6
local ccc = hc.curve(16384, 0, 
	{
		hc.segment(1/div, 1, hc.cubic_curve()),
		hc.segment(1/div, 0.8, hc.diocles_curve(1)),
		hc.segment(1/div, 1, hc.linear_curve()),
		hc.segment(1/div, 0.25, hc.cubic_curve()),
		hc.segment(1/div, 0.388, hc.cissoid_curve(1)),
		hc.segment(1/div, 0.123, hc.linear_curve())

	})

hc.write_as_wav("/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug/ttt.wav", ccc)
