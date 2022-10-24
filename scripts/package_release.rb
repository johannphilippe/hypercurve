require("fileutils")
require("os")
exec_extention = ""
extention = ""
prefix = ""
osname = ""
if(OS.linux?) then 
  osname = "linux"
  extension = ".so"
  prefix = "lib"
elsif(OS.windows?) then 
  osname = "windows"
  exec_extention = ".exe"
  extension = ".dll"
elsif(OS.mac?) then 
  osname = "macos"
  extension = ".dylib"
  prefix = "lib"
end


# Args : 1 path to build folder, 2 compiler name
puts "Preparing package for HYPERCURVE" 
path = ARGV[0]
compiler_name = ARGV[1]


# Platform name 
package_name = "HYPERCURVE_" + osname + "_" + compiler_name
outfolder = File.join(path, package_name)

puts "input path is " + path
puts "outfolder -> " +  outfolder

FileUtils.mkdir_p(outfolder)

Dir.foreach(path) do |filename|
  next if filename == '.' or filename == '..'

  puts filename
  if(filename == "faust_lib") then 
    FileUtils.cp_r( File.join(path, filename), outfolder)
  elsif(filename == "bin") then 
    binpath = File.join(path, "bin")
    Dir.foreach(binpath) do |libname|
      next if libname == '.' or libname == '..'
      if(libname =~ /csound_hypercurve#{extension}/) then 
        cs_path = File.join(outfolder, "Csound")
        FileUtils.mkdir_p(cs_path)
        FileUtils.cp( File.join(binpath, libname), cs_path)
      elsif (libname =~ /lua_hypercurve#{extension}/) then
        lua_path = File.join(outfolder, "Lua")
        FileUtils.mkdir_p(lua_path)
        FileUtils.cp( File.join(binpath, libname), lua_path)
      elsif (libname =~ /luajit#{extension}/) then 
        lua_path = File.join(outfolder, "Lua")
        FileUtils.mkdir_p(lua_path)
        FileUtils.cp( File.join(binpath, libname), lua_path)
      elsif (libname == "#{prefix}hypercurve#{extension}") then
        puts "hypercurve shared"
        hc_path = File.join(outfolder, "Hypercurve")
        FileUtils.mkdir_p(hc_path)
        FileUtils.cp( File.join(binpath, libname), hc_path)
      elsif (libname =~ /luajit#{exec_extention}/) then
        lua_path = File.join(outfolder, "Lua")
        FileUtils.mkdir_p(lua_path)
        FileUtils.cp( File.join(binpath, libname), lua_path)
      elsif(libname =~ /sndfile/) then
        lua_path = File.join(outfolder, "Lua")
        hc_path = File.join(outfolder, "Hypercurve")
        FileUtils.mkdir_p(lua_path)
        FileUtils.mkdir_p(hc_path)
        FileUtils.cp( File.join(binpath, libname), lua_path)
        FileUtils.cp( File.join(binpath, libname), hc_path)
      end
    end
  end
end

