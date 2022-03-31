#include"lauxlib.h"
#include"lua.h"
#include<unistd.h>
#include<memory.h>

#include"../src/core.h"
#include"../src/curve_lib.h"

// From lua 5.2
static void luaL_setfuncs (lua_State *L, const luaL_Reg *l, int nup) {
  luaL_checkstack(L, nup+1, "too many upvalues");
  for (; l->name != NULL; l++) {  /* fill the table with given functions */
    int i;
    lua_pushstring(L, l->name);
    for (i = 0; i < nup; i++)  /* copy upvalues to the top */
      lua_pushvalue(L, -(nup+1));
    lua_pushcclosure(L, l->func, nup);  /* closure with those upvalues */
    lua_settable(L, -(nup + 3));
  }
  lua_pop(L, nup);  /* remove upvalues */
}


#define luaL_newlibtable(L,l)   \
  lua_createtable(L, 0, sizeof(l)/sizeof((l)[0]) - 1)

#define luaL_newlib(L,l)        (luaL_newlibtable(L,l), luaL_setfuncs(L,l,0))

struct luahc_curve_base_t
{
    luahc_curve_base_t(hypercurve::curve_base *cb)
        : crv( cb )
    {}
    std::shared_ptr<hypercurve::curve_base> crv;
};

struct luahc_bezier_curve_base_t
{
    luahc_bezier_curve_base_t(hypercurve::bezier_curve_base *cb)
        : crv( cb )
    {}
    std::shared_ptr<hypercurve::bezier_curve_base> crv;
};

static luahc_curve_base_t** lua_curve_helper(lua_State *lua, hypercurve::curve_base *cb)
{
    luahc_curve_base_t **crv = (luahc_curve_base_t **) lua_newuserdata(lua, sizeof(luahc_curve_base_t*));
    *crv = new luahc_curve_base_t(cb);
    return crv;
}

static luahc_bezier_curve_base_t** lua_bezier_curve_helper(lua_State *lua, hypercurve::bezier_curve_base *cb)
{
    luahc_bezier_curve_base_t **crv = (luahc_bezier_curve_base_t **) lua_newuserdata(lua, sizeof(luahc_bezier_curve_base_t*));
    *crv = new luahc_bezier_curve_base_t(cb);
    return crv;
}

//extern "C" {
///////////////////////////////////////////:
// Control point
///////////////////////////////////////////:
int luahc_control_point(lua_State *lua)
{
    double x = lua_tonumber(lua, 1);
    double y = lua_tonumber(lua, 2);

    hypercurve::control_point *ctp = new hypercurve::control_point(x, y);
    hypercurve::control_point **cp = (hypercurve::control_point **)lua_newuserdata(lua, sizeof(hypercurve::control_point*));
    *cp = ctp;
    luaL_getmetatable(lua, "hypercurve.control_point");
    lua_setmetatable(lua, -2);
    return 1;
}

int luahc_control_point_x(lua_State *lua)
{
    int args = lua_gettop(lua);
    hypercurve::control_point *ctp = *(hypercurve::control_point **) luaL_checkudata(lua, 1, "hypercurve.control_point");
    if(args > 1)
    {
        double x_ = lua_tonumber(lua, 2);
        ctp->x = x_;
    }
    lua_pushnumber(lua, ctp->x);
    return 1;
}

int luahc_control_point_y(lua_State *lua)
{
    int args = lua_gettop(lua);
    hypercurve::control_point *ctp = *(hypercurve::control_point **) luaL_checkudata(lua, 1, "hypercurve.control_point");
    if(args > 1)
    {
        double y_ = lua_tonumber(lua, 2);
        ctp->y = y_;
    }
    lua_pushnumber(lua, ctp->y);
    return 1;
}

int luahc_control_point_xy(lua_State *lua)
{
    int args = lua_gettop(lua);
    hypercurve::control_point *ctp = *(hypercurve::control_point **) luaL_checkudata(lua, 1, "hypercurve.control_point");
    if(args > 2)
    {
        double x_ = lua_tonumber(lua, 2);
        double y_  = lua_tonumber(lua, 3);
        ctp->x = x_;
        ctp->y = y_;
    }
    lua_pushnumber(lua, ctp->x);
    lua_pushnumber(lua, ctp->y);
    return 2;
}

///////////////////////////////////////////:
// Create curve element
// Syntax : hypercurve.curve_base()
///////////////////////////////////////////:
int luahc_curve_base(lua_State *lua)
{
    luahc_curve_base_t **crv = lua_curve_helper(lua, new hypercurve::curve_base);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

int luahc_cubic_curve(lua_State *lua)
{
    luahc_curve_base_t **crv = lua_curve_helper(lua, new hypercurve::cubic_curve);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

int luahc_cissoid_curve(lua_State *lua)
{
    const double a = lua_tonumber(lua, 1);
    luahc_curve_base_t **crv = lua_curve_helper(lua, new hypercurve::cissoid_curve(a));
    //hypercurve::cissoid_curve **crv = (hypercurve::cissoid_curve **) lua_newuserdata(lua, sizeof(hypercurve::cissoid_curve*));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

int luahc_quadratic_bezier_curve(lua_State *lua)
{
    hypercurve::control_point *cp = *(hypercurve::control_point**) luaL_checkudata(lua, 1, "hypercurve.control_point");
    luahc_bezier_curve_base_t **crv = lua_bezier_curve_helper(lua, new hypercurve::quadratic_bezier_curve(*cp));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

int luahc_cubic_bezier_curve(lua_State *lua)
{
    hypercurve::control_point *cp = *(hypercurve::control_point**) luaL_checkudata(lua, 1, "hypercurve.control_point");
    hypercurve::control_point *cp2 = *(hypercurve::control_point**) luaL_checkudata(lua, 2, "hypercurve.control_point");
    luahc_bezier_curve_base_t **crv = lua_bezier_curve_helper(lua, new hypercurve::cubic_bezier_curve(*cp, *cp2));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

///////////////////////////////////////////:
// Create curve segment
// Syntax : hypercurve.segment(double frac, double y_dest, curve crv)
///////////////////////////////////////////:
int luahc_segment(lua_State *lua)
{
    double frac = lua_tonumber(lua, 1);
    double y_dest = lua_tonumber(lua, 2);
    luahc_curve_base_t *curve = *(luahc_curve_base_t **)luaL_checkudata(lua, 3, "hypercurve.curve_base");
    //luahc_curve_base_t *curve = *(luahc_curve_base_t **)lua_touserdata(lua, 3);
    // Memory --> TO FIX

    std::shared_ptr<hypercurve::segment> **l_seg = (std::shared_ptr<hypercurve::segment> **)
            lua_newuserdata(lua, sizeof(std::shared_ptr<hypercurve::segment>*));
    *l_seg = new std::shared_ptr<hypercurve::segment>(std::make_shared<hypercurve::segment>(frac, y_dest, curve->crv));

    luaL_getmetatable(lua, "hypercurve.segment");
    lua_setmetatable(lua, -2);

    return 1;
}

///////////////////////////////////////////:
// Create curve
// Syntax : hypercurve.curve()
///////////////////////////////////////////:

int luahc_curve(lua_State *lua)
{
    int definition = lua_tointeger(lua, 1);
    double y_start = lua_tonumber(lua, 2);
    luaL_checktype(lua, 3, LUA_TTABLE);
    size_t size = lua_objlen(lua, 3);

    std::vector< std::shared_ptr<hypercurve::segment> > segs;
    std::shared_ptr<hypercurve::segment> *seg_ptr = segs.data();
    for(size_t i = 1; i <= size; i++)
    {
        lua_rawgeti(lua, 3, i);
        std::shared_ptr<hypercurve::segment> *seg = *(std::shared_ptr<hypercurve::segment> **)luaL_checkudata(lua, -1, "hypercurve.segment");
        segs.push_back(*seg);
        seg_ptr++;
        lua_pop(lua, 1);
    }

    hypercurve::curve *crv = new hypercurve::curve(definition, y_start, segs);
    hypercurve::curve **l_crv = (hypercurve::curve **) lua_newuserdata(lua, sizeof(hypercurve::curve *));
    *l_crv = crv;
    luaL_getmetatable(lua, "hypercurve.curve");
    lua_setmetatable(lua, -2);

    return 1;
}


///////////////////////////////////////////:
// Curve to Soundfile
// Syntax : hypercurve.write_as_wav(string path, curve c)
///////////////////////////////////////////:
#include"sndfile.hh"
int luahc_curve_write_as_wav(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    std::string path = lua_tostring(lua, 2);
    SndfileHandle sf(path , SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_PCM_24, 1, 48000);

    sf.writef(crv->get_samples(), crv->get_definition());
    return 0;
}

int luahc_curve_ascii_display(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    std::string name = lua_tostring(lua, 2);
    std::string label = lua_tostring(lua, 3);
    std::string c = lua_tostring(lua, 4);
    crv->ascii_display(name, label, c[0]);
    return 0;
}

extern "C" {
///////////////////////////////////////////:
// Methods registration
///////////////////////////////////////////:

static int luahc_curve_class_gc(lua_State *lua) {
    //printf("## __gc\n");
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    delete crv;
    return 0;
}

static int luahc_index(lua_State *L) {
    //printf("## index\n");
    int i = luaL_checkinteger(L, 2);
    lua_pushinteger(L, i);
    return 1;
}

static int luahc_newindex(lua_State *L) {
    //printf("## index\n");
    return 0;
}

static const luaL_Reg luahc_curve_class_meta[] =
{
    { "__gc"        ,luahc_curve_class_gc          },
    //{ "__newindex"        ,luahc_newindex          },
    //{ "__index"        ,luahc_index          },
    { NULL          ,NULL            }
};
static const luaL_Reg luahc_curve_class_meth[] =
{
    { "ascii_display" ,luahc_curve_ascii_display },
    { "write_as_wav" ,luahc_curve_write_as_wav },
    { NULL          ,NULL            }
};

static const luaL_Reg luahc_static_meta[] =
{
    { "__index" ,luahc_index },
    //{ "__gc"        ,luahc_curve_class_gc          },
    //{ "__newindex"        ,luahc_newindex          },
    //{ "__call"  ,luahc_curve   },
    { NULL      ,NULL      }
};
static const luaL_Reg luahc_static_meth[] =
{
    {"new" , luahc_curve },
    {"segment", luahc_segment},
    {"control_point", luahc_control_point},

    {"curve_base", luahc_curve_base},
    {"linear", luahc_curve_base},
    {"cubic", luahc_cubic_curve},
    {"cissoid", luahc_cissoid_curve},
    {"diocles", luahc_cissoid_curve},
    {"quadratic_bezier", luahc_quadratic_bezier_curve},
    {"cubic_bezier", luahc_cubic_bezier_curve},
    { NULL      ,NULL      }
};


int luaopen_liblua_hypercurve (lua_State *lua)
{
    luaL_newmetatable(lua, "hypercurve.segment");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);
    luaL_newmetatable(lua, "hypercurve.curve_base");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);
    luaL_newmetatable(lua, "hypercurve.control_point");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);

    luaL_newmetatable(lua, "hypercurve.curve");
    luaL_setfuncs(lua, luahc_curve_class_meta, 0);
    luaL_newlib(lua, luahc_curve_class_meth);
    lua_setfield(lua, -2, "__index");
    lua_pop(lua, 1);

    luaL_newlib      (lua, luahc_static_meth);
    luaL_newlib      (lua, luahc_static_meta);
    lua_setmetatable (lua, -2);

    return 1;
}

}
