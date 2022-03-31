package.cpath = package.cpath .. ";/home/johann/Documents/GitHyb/build-hypercurve-Clang10-Debug"

local hc =  require("liblua_hypercurve")
local div = 6

function tp(v, str)

	for k,v in pairs(v) do
		print(str, k)
		if(type(v) == "table") then
			a = str .. "\t"
			tp(v, a)
		end
	end
end
print("Print Lua module components")
tp(hc, "")


-- Here a simple composite curve example
local ccc = hc.new(16384, 0, 
	{
		hc.segment(1/div, 1, hc.cubic()),
		hc.segment(1/div, 0.8, hc.diocles(1)),
		hc.segment(1/div, 1, hc.linear()),
		hc.segment(1/div, 0.25, hc.cubic()),
		hc.segment(1/div, 0.388, hc.cissoid(1)),
		hc.segment(1/div, 0.123, hc.linear())

	})
	

	ccc:ascii_display("luacurve", "y = luacurve(x)", "*")


--ccc:write_as_wav("/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug/ttt.wav")


-- Then a bezier curve example
local bez = hc.new(16384, 0, 
	{
		hc.segment(1/2, 1, hc.quadratic_bezier( 
					hc.control_point(1/4, 0.1))),
		hc.segment(1/2, 0, hc.cubic_bezier(
					hc.control_point(1/4, 0.9), hc.control_point(0.5, 0.9))),
	})

--hc.write_as_wav("/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug/bezier_lua.wav", bez)
bez:ascii_display("bezier", "y=bezier", "-")
