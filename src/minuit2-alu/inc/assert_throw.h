//
// define convenient macro to replace assertions with exceptions
//

#ifndef assert_throw_h
#define assert_throw_h

#define S(x) #x
#define S_(x) S(x)
#define S__LINE__ S_(__LINE__)
#define assert_throw(expr, msg) \
   (!(expr) \
    ? __ASSERT_VOID_CAST (0) \
    : throw std::runtime_error( __FILE__ "("  S__LINE__ "): " msg))

#endif
