-- Path to liblua_hypercurve
package.cpath = package.cpath .. ";/home/johann/Documents/GitHub/build-hypercurve-Clang10-Debug;"

local hc =  require("liblua_hypercurve")

-- This script will create a curve and write it in automation for every selected item in current Reaper project

local item_count = reaper.CountSelectedMediaItems(0)
if(item_count <= 0) then return end

local def = 1024
local crv = hc.curve(def, 0, hc.segment(1, 1, hc.diocles(1)))

for i = 0, item_count - 1 do
	local item = reaper.GetSelectedMediaItem(0, i)
	local duration = reaper.GetMediaItemInfo_Value(item, "D_LENGTH")
	local pos = reaper.GetMediaItemInfo_Value(item, "D_POSITION")
	local grain = duration / def
	local track = reaper.GetMediaItemTrack( item )
	local env = reaper.GetTrackEnvelopeByName(track, "Volume")
	
	for j = 0, def - 1 do
		local sample = crv:get_sample_at(j)
		reaper.InsertEnvelopePoint(env, pos + (grain * j), sample, 0, 0, false, false)
	end

 --reaper.SetEnvelopePoint( envelope, ptidx, timeIn, valueIn, shapeIn, tensionIn, selectedIn, noSortIn )
 -- reaper.InsertEnvelopePoint( envelope, time, value, shape, tension, selected, noSortIn )
 --  reaper.GetTrackEnvelopeByName( track, envname )
end

