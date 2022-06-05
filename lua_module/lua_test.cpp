// doscript.c
#include<iostream>

extern "C" {
#include "lauxlib.h"
#include "lua.h"
#include "lualib.h"

void isError(int error) {
    if(error) std::cout << "Lua error in script. ERROR CODE = " << error << std::endl;
       if(error != 0) {
               switch(error) {
               case LUA_ERRRUN:
                   std::cout << "Lua runtime error" << std::endl;
                   break;
               case LUA_ERRERR:
                   std::cout << "Lua handling function error " << std::endl;
                   break;
               case LUA_ERRMEM:
                   std::cout << "Lua memory allocation error " << std::endl;
                   break;
               case LUA_ERRFILE:
                   std::cout << "Lua Error File " << std::endl;
                   break;
               case LUA_ERRSYNTAX:
                   std::cout << "Lua syntax error " << std::endl;
                   break;
               case LUA_YIELD:
                   std::cout << "Lua yield error " << std::endl;
                   break;
               }
       }
}

static void fatal(const char* message) {
  fprintf(stderr,"%s\n", message);
  exit(EXIT_FAILURE);
}

#define LUA_OK 0

int main() {
  std::cout << "Create lua state " << std::endl;
  lua_State *L = luaL_newstate();
  std::cout << "opening libs " << std::endl;
  luaL_openlibs(L);
  std::cout << "calling test.lua" << std::endl;
  int res = luaL_dofile(L, "test.lua");
  if (res != LUA_OK) {
    fatal(lua_tostring(L,-1));
  }
  std::cout << "result = " << res << std::endl;
  isError(res);
  return res;
}

}
