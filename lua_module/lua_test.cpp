// doscript.c
extern "C" {
#include "lauxlib.h"
#include "lua.h"
#include "lualib.h"
#ifdef _WIN32
#include"SHARED_EXPORTS.h"
#endif
}

int main() {
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);
  luaL_dofile(L, "test.lua");
  return 0;
}
