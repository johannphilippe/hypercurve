#include"lauxlib.h"
#include"lua.h"
#include<unistd.h>
#include<memory.h>

#include"../src/core.h"
#include"../src/curve_lib.h"

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
    luaL_getmetatable(lua, "hypercurve.bezier_curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

int luahc_cubic_bezier_curve(lua_State *lua)
{
    hypercurve::control_point *cp = *(hypercurve::control_point**) luaL_checkudata(lua, 1, "hypercurve.control_point");
    hypercurve::control_point *cp2 = *(hypercurve::control_point**) luaL_checkudata(lua, 2, "hypercurve.control_point");
    luahc_bezier_curve_base_t **crv = lua_bezier_curve_helper(lua, new hypercurve::cubic_bezier_curve(*cp, *cp2));
    luaL_getmetatable(lua, "hypercurve.bezier_curve_base");
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

int luahc_bezier_segment(lua_State *lua)
{
    double frac = lua_tonumber(lua, 1);
    double y_dest = lua_tonumber(lua, 2);
    luahc_bezier_curve_base_t *curve = *(luahc_bezier_curve_base_t **)luaL_checkudata(lua, 3, "hypercurve.bezier_curve_base");

    std::shared_ptr<hypercurve::bezier_segment> **l_seg = (std::shared_ptr<hypercurve::bezier_segment> **)
            lua_newuserdata(lua, sizeof(std::shared_ptr<hypercurve::bezier_segment>*));
    *l_seg = new std::shared_ptr<hypercurve::bezier_segment>(std::make_shared<hypercurve::bezier_segment>(frac, y_dest, curve->crv));

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
        std::cout << "i  " << i << std::endl;
        lua_rawgeti(lua, 3, i);
        std::shared_ptr<hypercurve::segment> *seg = *(std::shared_ptr<hypercurve::segment> **)luaL_checkudata(lua, -1, "hypercurve.segment");
        segs.push_back(*seg);
        seg_ptr++;
        lua_pop(lua, 1);
    }
    std::cout << "segs recuperes " << std::endl;

    hypercurve::curve *crv = new hypercurve::curve(definition, y_start, segs);
    std::cout << "crv created " << std::endl;
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
int luahc_write_as_wav(lua_State *lua)
{
    std::string path = lua_tostring(lua, 1);
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 2, "hypercurve.curve");
    SndfileHandle sf(path , SFM_WRITE, SF_FORMAT_WAV | SF_FORMAT_PCM_24, 1, 48000);

    sf.writef(crv->get_samples(), crv->get_definition());
    return 0;
}


extern "C" {
///////////////////////////////////////////:
// Methods registration
///////////////////////////////////////////:

static const luaL_Reg luahc_api[] = {
    // Basic getters and setters
    {"curve", luahc_curve},
    {"segment", luahc_segment},
    {"bezier_segment", luahc_bezier_segment},

    // Curves
    {"curve_base", luahc_curve_base},
    {"linear_curve", luahc_curve_base},
    {"cubic_curve", luahc_cubic_curve},
    {"cissoid_curve", luahc_cissoid_curve},
    {"diocles_curve", luahc_cissoid_curve},
    {"quadratic_bezier_curve", luahc_quadratic_bezier_curve},
    {"cubic_bezier_curve", luahc_cubic_bezier_curve},

    {"control_point", luahc_control_point},

    {"write_as_wav", luahc_write_as_wav},

    {NULL, NULL}
};

int luaopen_liblua_hypercurve (lua_State *lua)
{
    luaL_newmetatable(lua, "hypercurve.curve");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);
    luaL_newmetatable(lua, "hypercurve.segment");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);
    luaL_newmetatable(lua, "hypercurve.curve_base");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);
    luaL_newmetatable(lua, "hypercurve.bezier_curve_base");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);
    luaL_newmetatable(lua, "hypercurve.control_point");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);

    lua_newtable(lua);
    luaL_register(lua, "hypercurve", luahc_api);
    /* ... more push operations ... */
    return 1;
}

}
