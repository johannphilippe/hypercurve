# Lua Hypercurve

Lua Hypercurve is a frontend to hypercurve. 
It allows you to write hybrid curves from Lua language. 

To include the module to your environment just run the following (replace .so with .dylib or .dll)

```lua
package.cpath = package.cpath .. ";/your/path/to/hypercurve/?.so;"
local hc = require("liblua_hypercurve")
```

Then you can use Hypercurve.