#ifndef ASSERT_THROW_H
#define ASSERT_THROW_H

#include <stdexcept>
#include <cassert>

// 1. Always undefine first to ensure our logic takes precedence
#ifdef assert
  #undef assert
#endif

#ifdef NDEBUG
  // Release mode: define as a no-op that is optimized away
  // This behaves like the "normal" assert when NDEBUG is defined.
  #define assert(expr) static_cast<void>(0)
#else
  // 2. Stringification helpers with unique prefixes
  #define ASSERT_THROW_STR_HELPER_(x) #x
  #define ASSERT_THROW_STR_(x) ASSERT_THROW_STR_HELPER_(x)

  // 3. The R-friendly replacement (Debug/Check Mode)
  // Replaces "Minuit2 Error at" with "assert at" as requested.
  #define assert(expr) \
    ((expr) \
     ? static_cast<void>(0) \
     : throw std::runtime_error("assert at " __FILE__ ":" ASSERT_THROW_STR_(__LINE__) " - Assertion '" #expr "' failed."))
#endif

#endif
