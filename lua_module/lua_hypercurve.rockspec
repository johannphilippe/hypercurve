#!/usr/bin/env lua

package = 'lua-hypercurve'
version = '0.0-1'
source  = {
    url    = 'git://github.com/johannphilippe/hypercurve.git',
    branch = 'v1.0',
}
description = {
    summary  = "A Lua binding to Ryan Dahl's http request/response parser.",
    detailed = '',
    homepage = 'https://github.com/johannphilippe/hypercurve',
    license  = 'MIT', --as with Ryan's
}
dependencies = {
    'lua >= 5.1'
}
build    = {
    type = 'cmake',
    variables = {
        INSTALL_CMOD      = "$(LIBDIR)",
        CMAKE_BUILD_TYPE  = "$(CMAKE_BUILD_TYPE)",
        ["CFLAGS:STRING"] = "$(CFLAGS)",
        BUILD_LUA_MODUEL = "TRUE",
    },
}