-- Find package if not in standard Lua CPATH
package.cpath = package.cpath .. ";/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug/?.so;"

local hc =  require("liblua_hypercurve")



print("Print Lua module components")
function tp(v, str)

	for k,v in pairs(v) do
		print(str, k)
		if(type(v) == "table") then
			a = str .. "\t"
			tp(v, a)
		end
	end
end
tp(hc, "")

-- Here a simple hybrid curve example
local div = 6
local hybrid = hc.curve(16384, 0, 
	{
		hc.segment(1/div, 1, hc.cubic()),
		hc.segment(1/div, 0.8, hc.diocles(1)),
		hc.segment(1/div, 1, hc.linear()),
		hc.segment(1/div, 0.25, hc.cubic()),
		hc.segment(1/div, 0.388, hc.cissoid(1)),
		hc.segment(1/div, 0.123, hc.linear())

	})
	

-- You can access samples as a table
local samples = hybrid:get_samples()
print("samples retrieved : ", #samples)

-- Or access samples with getter
local smp = hybrid:get_sample_at(5)
print("single sample smp = ", smp) 

-- Display your curve in terminal
hybrid:ascii_display("luacurve", "y = luacurve(x)", "*")


--The following writes as wave file with specified number of points (definition)
--hybrid:write_as_wav("/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug/ttt.wav")

-- Then a bezier curve example
--
local bez = hc.curve(16384, 0, 
	{
		hc.segment(1/2, 1, hc.quadratic_bezier( 
					hc.control_point(1/4, 0.1))),
		hc.segment(1/2, 0, hc.cubic_bezier(
					hc.control_point(1/4, 0.9), hc.control_point(0.5, 0.9))),
	})

-- Control points also have class methods, that can be used as setters and/or getters
-- local cp = control_point(1,0.5)
-- local x = cp:x(1.5)
-- local y = cp:y(0.5)
-- local xx, yy = cp:xy(1.5, 2.5)

bez:ascii_display("bezier", "y=bezier", "-")
--hc:write_as_wav("/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug/bezier_lua.wav", bez)

--
-- User defined curve with user defined function

function cube(x)
	return x*x*x
end

local crv = hc.user_defined(cube)
print("ud ccc ok")
local seg = hc.segment(1, 1, crv)
print("ud seg ok ")
local full_curve = hc.curve(1024, 0, {seg} )
print("ud curve ok")
full_curve:ascii_display("user defined", "user(x)", "*")

-- other curve types that can be  passed as segment argument

-- Spline : takes a list of control points as argument (any number of control points)
local spl = hc.cubic_spline({hc.control_point(0.2, 0.6), hc.control_point(0.5, 0.1), hc.control_point(0.3, 0.345)})

-- Catmull Rom : takes two control points : P0 and P3 (P0.x must be negative, and P3.x must be more than 1)
local cm = hc.catmull_rom(hc.control_point(-2, 1), hc.control_point(3, 5))
local cm_seg = hc.segment(1, 1, cm)
local cm_crv = hc.curve( 16384, 0, 
	{
		cm_seg
	})

cm_crv:ascii_display("catmullrom", "cm(x)", "*")

