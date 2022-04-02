ARGV.each do |a|
	dir = File.dirname(a)
	name = File.basename(a, ".*")
	ext = File.extname(a)

	file_new = dir + "/" + name + ".png"
    
    cmd = "ffmpeg -i \"#{a}\" -filter_complex showwavespic=draw=full:filter=peak:s=1024x960 \"#{file_new}\""
    system(cmd)

    cmd = "convert \"#{file_new}\" -gravity north -crop 100x50% \"#{file_new}\""
    system(cmd)

end



