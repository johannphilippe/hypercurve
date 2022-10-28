/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

extern "C" {
#include"lauxlib.h"
#include"lua.h"
}

#ifndef WIN32
    #include<unistd.h>
    #define SHARED_EXPORTS
    #define LUA_EXPORT extern
#else
extern "C" {
    #include "SHARED_EXPORTS.h"
    #define LUA_EXPORT extern SHARED_EXPORTS
}
#endif

#include<memory.h>
#include<functional>

#include"../src/core.h"
#include"../src/curve_lib.h"
#include"../src/utilities.h"
#include"compat-5.3.h"

struct SHARED_EXPORTS luahc_curve_base_t
{
    luahc_curve_base_t(hypercurve::curve_base *cb)
        : crv( cb )
    {}
    std::shared_ptr<hypercurve::curve_base> crv;
};

LUA_EXPORT luahc_curve_base_t** lua_curve_helper(lua_State *lua, hypercurve::curve_base *cb)
{
    luahc_curve_base_t **crv = (luahc_curve_base_t **) lua_newuserdata(lua, sizeof(luahc_curve_base_t*));
    *crv = new luahc_curve_base_t(cb);
    return crv;
}

struct SHARED_EXPORTS luahc_user_defined_curve_t : public hypercurve::user_defined_curve
{
    luahc_user_defined_curve_t(int reg, lua_State *lua)
        : _lua_registry_index(reg)
        , lua_state(lua)
    {
        callback = [&](double x)
        {
            lua_rawgeti(lua_state, LUA_REGISTRYINDEX, _lua_registry_index );
            lua_pushnumber(lua_state, x);
            if ( 0 != lua_pcall( lua_state, 1, 1, 0 ) ) {
              printf("Failed to call the callback!\n %s\n", lua_tostring(lua_state, -1 ) );
              return 0.;
            }
            double res = lua_tonumber(lua_state, -1);
            return res;

        };
    }
    int _lua_registry_index;
    lua_State *lua_state;
};

///////////////////////////////////////////
// Control point
///////////////////////////////////////////
LUA_EXPORT int luahc_control_point(lua_State *lua)
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

LUA_EXPORT int luahc_control_point_x(lua_State *lua)
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

LUA_EXPORT int luahc_control_point_y(lua_State *lua)
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

LUA_EXPORT int luahc_control_point_xy(lua_State *lua)
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

///////////////////////////////////////////
// Create curve element
// Syntax : hypercurve.curve_base()
///////////////////////////////////////////

LUA_EXPORT int luahc_curve_base(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::curve_base);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_cubic_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::cubic_curve);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_logarithmic_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::logarithmic_curve);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_exponential_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::exponential_curve);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_power_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::power_curve(lua_tonumber(lua, 1)));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_cissoid_curve(lua_State *lua)
{
    const double a = lua_tonumber(lua, 1);
    lua_curve_helper(lua, new hypercurve::cissoid_curve(a));
    //hypercurve::cissoid_curve **crv = (hypercurve::cissoid_curve **) lua_newuserdata(lua, sizeof(hypercurve::cissoid_curve*));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_hanning_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::hanning_curve);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}
LUA_EXPORT int luahc_hamming_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::hamming_curve);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}
LUA_EXPORT int luahc_blackman_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::blackman_curve);
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_gauss_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::gauss_curve(
                         lua_tonumber(lua, 1), lua_tonumber(lua, 2)));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_catenary_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::catenary_curve(lua_tonumber(lua, 1)));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_toxoid_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::toxoid_curve(lua_tonumber(lua, 1)));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_mouse_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::mouth_curve());
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_bicorn_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::bicorn_curve(lua_toboolean(lua, 1)));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_typed_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::typed_curve(
                         lua_tonumber(lua, 1)));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_tightrope_walker_curve(lua_State *lua)
{
    lua_curve_helper(lua, new hypercurve::tightrope_walker_curve(
                         lua_tonumber(lua, 1), lua_tonumber(lua, 2)));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_quadratic_bezier_curve(lua_State *lua)
{
    hypercurve::control_point *cp = *(hypercurve::control_point**) luaL_checkudata(lua, 1, "hypercurve.control_point");
    lua_curve_helper(lua, new hypercurve::quadratic_bezier_curve(*cp));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_cubic_bezier_curve(lua_State *lua)
{
    hypercurve::control_point *cp = *(hypercurve::control_point**) luaL_checkudata(lua, 1, "hypercurve.control_point");
    hypercurve::control_point *cp2 = *(hypercurve::control_point**) luaL_checkudata(lua, 2, "hypercurve.control_point");
    lua_curve_helper(lua, new hypercurve::cubic_bezier_curve(*cp, *cp2));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_cubic_spline_curve(lua_State *lua)
{
    luaL_checktype(lua, 1, LUA_TTABLE);
    size_t size = luaL_len(lua, 1); //lua_objlen(lua, 3);
    std::vector<hypercurve::control_point> cps;
    for(size_t i = 1; i <= size; i++)
    {
        lua_rawgeti(lua, 1, i);
        hypercurve::control_point *cp = *(hypercurve::control_point **) luaL_checkudata(lua, -1, "hypercurve.control_point");
        cps.push_back(*cp);
        lua_pop(lua, 1);
    }
    lua_curve_helper(lua, new hypercurve::cubic_spline_curve(cps));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_catmull_rom_spline_curve(lua_State *lua)
{
    hypercurve::control_point *cp = *(hypercurve::control_point**) luaL_checkudata(lua, 1, "hypercurve.control_point");
    hypercurve::control_point *cp2 = *(hypercurve::control_point**) luaL_checkudata(lua, 2, "hypercurve.control_point");
    lua_curve_helper(lua, new hypercurve::catmull_rom_spline_curve(*cp, *cp2));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_user_defined_curve(lua_State *lua)
{
    if(lua_gettop(lua) == 1 && lua_isfunction(lua, -1))
    {
        int callback_reference = luaL_ref( lua, LUA_REGISTRYINDEX );
        lua_curve_helper(lua, new luahc_user_defined_curve_t(callback_reference, lua));
        luaL_getmetatable(lua, "hypercurve.curve_base");
        lua_setmetatable(lua, -2);
        return 1;
    }
    return 0;
}

LUA_EXPORT int luahc_lagrange_interpolation_curve(lua_State *lua)
{
    hypercurve::memory_vector< hypercurve::control_point> vec(lua_gettop(lua));
    for(size_t i = 1; i <= vec.size(); ++i)
    {
        hypercurve::control_point *cp = *(hypercurve::control_point**) luaL_checkudata(lua, i, "hypercurve.control_point");
        vec[i - 1] = *cp;
    }
    lua_curve_helper(lua, new hypercurve::lagrange_polynomial_curve(vec));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_polynomial_curve(lua_State *lua)
{
    hypercurve::memory_vector<double> vec(lua_gettop(lua));
    for(size_t i = 1; i <= vec.size(); ++i)
        vec[i-1] = lua_tonumber(lua, i);
    lua_curve_helper(lua, new hypercurve::polynomial_curve(vec));
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

///////////////////////////////////////////
// Invert curve
///////////////////////////////////////////

LUA_EXPORT int luahc_invert_curve(lua_State *lua)
{
    luahc_curve_base_t *curve = *(luahc_curve_base_t **)luaL_checkudata(lua, 1, "hypercurve.curve_base");
    curve->crv->inverted = true;
    luahc_curve_base_t **crv = (luahc_curve_base_t **) lua_newuserdata(lua, sizeof(luahc_curve_base_t*));
    *crv = curve;
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_mirror_curve(lua_State *lua)
{
    luahc_curve_base_t *curve = *(luahc_curve_base_t **)luaL_checkudata(lua, 1, "hypercurve.curve_base");
    curve->crv->mirrored = true;
    luahc_curve_base_t **crv = (luahc_curve_base_t **) lua_newuserdata(lua, sizeof(luahc_curve_base_t*));
    *crv = curve;
    luaL_getmetatable(lua, "hypercurve.curve_base");
    lua_setmetatable(lua, -2);
    return 1;
}

///////////////////////////////////////////
// Create curve segment
// Syntax : hypercurve.segment(double frac, double y_dest, curve crv)
///////////////////////////////////////////

LUA_EXPORT int luahc_segment(lua_State *lua)
{
    double frac = lua_tonumber(lua, 1);
    double y_dest = lua_tonumber(lua, 2);
    luahc_curve_base_t *curve = *(luahc_curve_base_t **)luaL_checkudata(lua, 3, "hypercurve.curve_base");
    //luahc_curve_base_t *curve = *(luahc_curve_base_t **)lua_touserdata(lua, 3);
    // Memory --> TO FIX ?

    hypercurve::segment **l_seg = (hypercurve::segment **)
            lua_newuserdata(lua, sizeof(hypercurve::segment*));
    *l_seg = new hypercurve::segment(frac, y_dest, curve->crv);

    luaL_getmetatable(lua, "hypercurve.segment");
    lua_setmetatable(lua, -2);

    return 1;
}

///////////////////////////////////////////:
// Create curve
// Syntax : hypercurve.curve()
///////////////////////////////////////////:

LUA_EXPORT int luahc_curve(lua_State *lua)
{
    int definition = lua_tointeger(lua, 1);
    double y_start = lua_tonumber(lua, 2);
    luaL_checktype(lua, 3, LUA_TTABLE);
    size_t size = luaL_len(lua, 3); //lua_objlen(lua, 3);

    std::vector< hypercurve::segment > segs;
    hypercurve::segment *seg_ptr = segs.data();
    for(size_t i = 1; i <= size; i++)
    {
        lua_rawgeti(lua, 3, i);
        hypercurve::segment *seg = *(hypercurve::segment **)luaL_checkudata(lua, -1, "hypercurve.segment");
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


LUA_EXPORT int luahc_concatenate(lua_State *lua)
{

    luaL_checktype(lua, 2, LUA_TTABLE);
    size_t definition = lua_tointeger(lua, 1);
    size_t arr_size = luaL_len(lua, 2);
    std::vector<hypercurve::curve *> to_concat;
    for(size_t i = 1; i <= arr_size; i++)
    {
        lua_rawgeti(lua, 2, i);
        hypercurve::curve *crv = *(hypercurve::curve **)luaL_checkudata(lua, -1, "hypercurve.curve");
        to_concat.push_back(crv);
        lua_pop(lua, 1);
    }

    hypercurve::curve *crv = new hypercurve::curve(definition, to_concat);
    hypercurve::curve **l_crv = (hypercurve::curve **)lua_newuserdata(lua, sizeof(hypercurve::curve *));
    *l_crv = crv;
    luaL_getmetatable(lua, "hypercurve.curve");
    lua_setmetatable(lua, -2);
    return 1;
}


///////////////////////////////////////////:
// Curve methods
///////////////////////////////////////////:
LUA_EXPORT int luahc_curve_ascii_display(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    std::string name = lua_tostring(lua, 2);
    std::string label = lua_tostring(lua, 3);
    std::string c = lua_tostring(lua, 4);
    crv->ascii_display(name, label, c[0]);
    return 0;
}

LUA_EXPORT int luahc_curve_write_png(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    std::string path = lua_tostring(lua, 2);
    int args = lua_gettop(lua);

    bool waveform = (args > 2) ? lua_toboolean(lua, 3) : false;
    bool fill = (args > 3) ? lua_toboolean(lua, 4) : true;
    bool draw_grid = (args  >  4) ? lua_toboolean(lua, 5) : true;
    bool invert = (args > 5) ? lua_toboolean(lua, 6) : false;
    hypercurve::png p(2048, 1024, invert ? hypercurve::white : hypercurve::black, invert ? hypercurve::red : hypercurve::purple);
    p.draw_curve(crv->get_samples(), crv->get_definition(), fill, waveform);
    if(draw_grid) p.draw_grid(10, 10, invert ? hypercurve::black : hypercurve::white);
    p.write_as_png(path);
    return 0;
}

// Return one sample
LUA_EXPORT int luahc_curve_get_sample_at(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    int index = lua_tonumber(lua, 2);
    lua_pushnumber(lua, crv->get_samples()[index]);
    return 1;
}

// Returns a lua table
LUA_EXPORT int luahc_curve_get_samples(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    lua_newtable(lua);
    for(size_t i = 0; i < crv->get_definition(); ++i)
    {
        lua_pushnumber(lua, crv->get_sample_at(i));
        lua_rawseti(lua, -2, i+1);
    }
    return 1;
}

LUA_EXPORT int luahc_curve_normalize_y(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    crv->normalize_y(lua_tonumber(lua, 2), lua_tonumber(lua, 3));
    return 0;
}

//extern "C" {
///////////////////////////////////////////:
// Methods registration
///////////////////////////////////////////:

LUA_EXPORT int luahc_curve_class_gc(lua_State *lua) {
    //printf("## __gc\n");
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    delete crv;
    return 0;
}

LUA_EXPORT int luahc_curve_class_add(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    hypercurve::curve *crv2 = *(hypercurve::curve **) luaL_checkudata(lua, 2, "hypercurve.curve");
    hypercurve::curve *res = new hypercurve::curve((*crv) + (*crv2));
    hypercurve::curve **l_crv = (hypercurve::curve **) lua_newuserdata(lua, sizeof(hypercurve::curve *));
    *l_crv = res;
    luaL_getmetatable(lua, "hypercurve.curve");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_curve_class_sub(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    hypercurve::curve *crv2 = *(hypercurve::curve **) luaL_checkudata(lua, 2, "hypercurve.curve");
    hypercurve::curve *res = new hypercurve::curve((*crv) - (*crv2));
    hypercurve::curve **l_crv = (hypercurve::curve **) lua_newuserdata(lua, sizeof(hypercurve::curve *));
    *l_crv = res;
    luaL_getmetatable(lua, "hypercurve.curve");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_curve_class_mult(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    hypercurve::curve *crv2 = *(hypercurve::curve **) luaL_checkudata(lua, 2, "hypercurve.curve");
    hypercurve::curve *res = new hypercurve::curve((*crv) * (*crv2));
    hypercurve::curve **l_crv = (hypercurve::curve **) lua_newuserdata(lua, sizeof(hypercurve::curve *));
    *l_crv = res;
    luaL_getmetatable(lua, "hypercurve.curve");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_curve_class_div(lua_State *lua)
{
    hypercurve::curve *crv = *(hypercurve::curve **) luaL_checkudata(lua, 1, "hypercurve.curve");
    hypercurve::curve *crv2 = *(hypercurve::curve **) luaL_checkudata(lua, 2, "hypercurve.curve");
    hypercurve::curve *res = new hypercurve::curve((*crv) / (*crv2));
    hypercurve::curve **l_crv = (hypercurve::curve **) lua_newuserdata(lua, sizeof(hypercurve::curve *));
    *l_crv = res;
    luaL_getmetatable(lua, "hypercurve.curve");
    lua_setmetatable(lua, -2);
    return 1;
}

LUA_EXPORT int luahc_control_point_class_gc(lua_State *lua)
{
    hypercurve::control_point *cp = *(hypercurve::control_point **) luaL_checkudata(lua, 1, "hypercurve.control_point");
    delete cp;
    return 0;
}

LUA_EXPORT int luahc_index(lua_State *L) {
    //printf("## index\n");
    int i = luaL_checkinteger(L, 2);
    lua_pushinteger(L, i);
    return 1;
}

LUA_EXPORT const luaL_Reg luahc_curve_class_meta[] =
{
    { "__gc"        ,luahc_curve_class_gc          },
    { "__add" 		,luahc_curve_class_add		   },
    { "__sub" 		,luahc_curve_class_sub		   },
    { "__mul" 		,luahc_curve_class_mult		   },
    { "__div" 		,luahc_curve_class_div		   },
    { NULL          ,NULL            }
};

LUA_EXPORT const luaL_Reg luahc_control_point_class_meta[] =
{
    { "__gc"        ,luahc_control_point_class_gc          },
    { NULL          ,NULL            }
};

LUA_EXPORT const luaL_Reg luahc_curve_class_meth[] =
{
    {"ascii_display" ,luahc_curve_ascii_display },
    {"write_as_png" ,luahc_curve_write_png },
    {"normalize", luahc_curve_normalize_y},
    {"get_samples", luahc_curve_get_samples},
    {"get_sample_at", luahc_curve_get_sample_at},
    // deprecated
    {"normalize_y", luahc_curve_normalize_y},
    { NULL          ,NULL            }
};

LUA_EXPORT const luaL_Reg luahc_control_point_class_meth[] =
{
    {"x", luahc_control_point_x},
    {"y", luahc_control_point_y},
    {"xy", luahc_control_point_xy},
    { NULL          ,NULL            }
};

LUA_EXPORT const luaL_Reg luahc_static_meta[] =
{
    { "__index" ,luahc_index },
    //{ "__gc"        ,luahc_curve_class_gc          },
    //{ "__newindex"        ,luahc_newindex          },
    //{ "__call"  ,luahc_curve   },
    { NULL      ,NULL      }
};

LUA_EXPORT const luaL_Reg luahc_static_meth[] =
{
    {"hypercurve" , luahc_curve },
    {"curve" , luahc_curve },  // alias for hypercurve

    {"segment", luahc_segment},

    {"control_point", luahc_control_point},
    {"point", luahc_control_point},

    {"curve_base", luahc_curve_base},
    {"linear", luahc_curve_base},
    {"cubic", luahc_cubic_curve},
    {"power", luahc_power_curve},
    {"diocles", luahc_cissoid_curve},
    {"hanning", luahc_hanning_curve},
    {"hamming", luahc_hamming_curve},
    {"blackman", luahc_blackman_curve},
    {"gauss", luahc_gauss_curve},
    {"toxoid", luahc_toxoid_curve},
    {"mouse", luahc_mouse_curve},
    {"bicorn", luahc_bicorn_curve},
    {"catenary", luahc_catenary_curve},
    {"tightrope_walker", luahc_tightrope_walker_curve},
    {"quadratic_bezier", luahc_quadratic_bezier_curve},
    {"cubic_bezier", luahc_cubic_bezier_curve},
    {"cubic_spline", luahc_cubic_spline_curve},
    {"catmull_rom", luahc_catmull_rom_spline_curve},
    {"lagrange_polynomial", luahc_lagrange_interpolation_curve},
    {"polynomial", luahc_polynomial_curve},
    {"typed", luahc_typed_curve},
    {"user_defined", luahc_user_defined_curve},

    // Add curve suffix
    {"linear_curve", luahc_curve_base},
    {"cubic_curve", luahc_cubic_curve},
    {"logarithmic_curve", luahc_logarithmic_curve},
    {"exponential_curve", luahc_exponential_curve},
    {"power_curve", luahc_power_curve},
    {"diocles_curve", luahc_cissoid_curve},
    {"hanning_curve", luahc_hanning_curve},
    {"hamming_curve", luahc_hamming_curve},
    {"blackman_curve", luahc_blackman_curve},
    {"gauss_curve", luahc_gauss_curve},
    {"toxoid_curve", luahc_toxoid_curve},
    {"mouth_curve", luahc_mouse_curve},
    {"bicorn_curve", luahc_bicorn_curve},
    {"catenary_curve", luahc_catenary_curve},
    {"tightrope_walker_curve", luahc_tightrope_walker_curve},
    {"quadratic_bezier_curve", luahc_quadratic_bezier_curve},
    {"cubic_bezier_curve", luahc_cubic_bezier_curve},
    {"cubic_spline_curve", luahc_cubic_spline_curve},
    {"catmull_rom_curve", luahc_catmull_rom_spline_curve},
    {"lagrange_polynomial_curve", luahc_lagrange_interpolation_curve},
    {"polynomial_curve", luahc_polynomial_curve},
    {"typed_curve", luahc_typed_curve},
    {"user_defined_curve", luahc_user_defined_curve},

    // Helpers
    {"invert", luahc_invert_curve},
    {"mirror", luahc_mirror_curve},
    {"concatenate", luahc_concatenate},

    // Aliases
    {"cissoid", luahc_cissoid_curve},
    {"gaussian", luahc_gauss_curve},
    {"duplicatrix_cubic", luahc_toxoid_curve},
    {"funicular", luahc_catenary_curve},
    {"kiss", luahc_mouse_curve},
    {"cocked_hat", luahc_bicorn_curve},

    // Add curve suffix for aliases
    {"cissoid_curve", luahc_cissoid_curve},
    {"gaussian_curve", luahc_gauss_curve},
    {"duplicatrix_cubic_curve", luahc_toxoid_curve},
    {"funicular_curve", luahc_catenary_curve},
    {"kiss_curve", luahc_mouse_curve},
    {"cocked_hat_curve", luahc_bicorn_curve},

    {"concat", luahc_concatenate},

    { NULL      ,NULL      }
};

extern "C" {
LUA_EXPORT int luaopen_lua_hypercurve(lua_State *lua)
{
    luaL_newmetatable(lua, "hypercurve.segment");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);
    luaL_newmetatable(lua, "hypercurve.curve_base");
    lua_pushstring(lua, "__index");
    lua_pushvalue(lua, -2);
    lua_settable(lua, -3);

    luaL_newmetatable(lua, "hypercurve.curve");
    luaL_setfuncs(lua, luahc_curve_class_meta, 0);
    luaL_newlib(lua, luahc_curve_class_meth);
    lua_setfield(lua, -2, "__index");
    lua_pop(lua, 1);

    luaL_newmetatable(lua, "hypercurve.control_point");
    luaL_setfuncs(lua, luahc_control_point_class_meta, 0);
    luaL_newlib(lua, luahc_control_point_class_meth);
    lua_setfield(lua, -2, "__index");
    lua_pop(lua, 1);

    luaL_newlib      (lua, luahc_static_meth);
    luaL_newlib      (lua, luahc_static_meta);
    lua_setmetatable (lua, -2);

    return 1;
}

}
