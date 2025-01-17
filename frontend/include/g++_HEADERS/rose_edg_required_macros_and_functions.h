/* JJW (12/8/2008): Do not include these definitions for EDG 3.10 */
/* #ifndef ROSE_USE_NEW_EDG_INTERFACE */
/* DQ (11/1/2011): I think we do need this for at least C++ code (e.g. __CHAR_BIT__ must
   be defined along with lots of other Common Predefined Macros.
*/
#if 1

/* DQ (6/12/2005): Permit this to be optionally commented out because it makes it hard to
    look at the AST merging graphs. The ROSE-IR graphs will now filter nodes associated
    with declarations from this file (marked as front-end specific).
*/
#if 0
#define SKIP_ROSE_BUILTIN_DECLARATIONS
#endif

/* This macro can be defined using "-DSKIP_ROSE_BUILTIN_DECLARATIONS" to avoid ROSE
   builtin functions required for compatability with the user selected backend compiler 
*/
#ifndef SKIP_ROSE_BUILTIN_DECLARATIONS

/* Must use C style comments so that "--edg:old_c" options will work! */
/* DQ (7/13/2006): Undefine these before defining them to avoid warnings. */
/* DQ (12/23/2006): Let EDG define this if possible, but we reset it for 64 bit systems
   where either EDG does not get it correct or we don't setup EDG correctly!
#undef __SIZE_TYPE__ 
*/
#undef __VERSION__

// $REPLACE_ME_WITH_MACRO_DEFINITIONS
#define __DBL_MIN_EXP__ (-1021)
#define __FLT_MIN__ 1.17549435e-38F
#define __CHAR_BIT__ 8
#define __WCHAR_MAX__ 2147483647
#define __GCC_HAVE_SYNC_COMPARE_AND_SWAP_1 1
#define __GCC_HAVE_SYNC_COMPARE_AND_SWAP_2 1
#define __GCC_HAVE_SYNC_COMPARE_AND_SWAP_4 1
#define __DBL_DENORM_MIN__ 4.9406564584124654e-324
#define __FLT_EVAL_METHOD__ 2
#define __unix__ 1
#define __DBL_MIN_10_EXP__ (-307)
#define __FINITE_MATH_ONLY__ 0
#define __DEC64_MAX_EXP__ 385
#define __SHRT_MAX__ 32767
#define __LDBL_MAX__ 1.18973149535723176502e+4932L
#define __UINTMAX_TYPE__ long long unsigned int
#define __linux 1
#define __DEC32_EPSILON__ 1E-6DF
#define __unix 1
#define __LDBL_MAX_EXP__ 16384
#define __linux__ 1
#define __SCHAR_MAX__ 127
#define __DBL_DIG__ 15
#define _FORTIFY_SOURCE 2
#define __SIZEOF_INT__ 4
#define __SIZEOF_POINTER__ 4
#define __USER_LABEL_PREFIX__ 
#define __LDBL_HAS_INFINITY__ 1
#define __FLT_EPSILON__ 1.19209290e-7F
#define __LDBL_MIN__ 3.36210314311209350626e-4932L
#define __DEC32_MAX__ 9.999999E96DF
#define __SIZEOF_LONG__ 4
#define __DECIMAL_DIG__ 21
#define __gnu_linux__ 1
#define __LDBL_HAS_QUIET_NAN__ 1
#define __GXX_RTTI 1
#define __FLT_HAS_DENORM__ 1
#define __SIZEOF_LONG_DOUBLE__ 12
#define __BIGGEST_ALIGNMENT__ 16
#define __DBL_MAX__ 1.7976931348623157e+308
#define __DBL_HAS_INFINITY__ 1
#define __DEC32_MIN_EXP__ (-94)
#define __LDBL_HAS_DENORM__ 1
#define __DEC128_MAX__ 9.999999999999999999999999999999999E6144DL
#define __DEC32_MIN__ 1E-95DF
#define __DEPRECATED 1
#define __DBL_MAX_EXP__ 1024
#define __DEC128_EPSILON__ 1E-33DL
#define __LONG_LONG_MAX__ 9223372036854775807LL
#define __SIZEOF_SIZE_T__ 4
#define __SIZEOF_WINT_T__ 4
#define __GCC_HAVE_DWARF2_CFI_ASM 1
#define __GXX_ABI_VERSION 1002
#define __FLT_MIN_EXP__ (-125)
#define __DBL_MIN__ 2.2250738585072014e-308
#define __FLT_MIN_10_EXP__ (-37)
#define __DECIMAL_BID_FORMAT__ 1
#define __DEC128_MIN__ 1E-6143DL
#define __REGISTER_PREFIX__ 
#define __DBL_HAS_DENORM__ 1
#define __NO_INLINE__ 1
#define __i386 1
#define __FLT_MANT_DIG__ 24
#define __VERSION__ "4.4.3"
#define __DEC64_EPSILON__ 1E-15DD
#define __DEC128_MIN_EXP__ (-6142)
#define __i486__ 1
#define unix 1
#define __i386__ 1
#define __SIZE_TYPE__ unsigned int
#define __ELF__ 1
#define __FLT_RADIX__ 2
#define __LDBL_EPSILON__ 1.08420217248550443401e-19L
#define __SIZEOF_PTRDIFF_T__ 4
#define __DEC32_SUBNORMAL_MIN__ 0.000001E-95DF
#define __FLT_HAS_QUIET_NAN__ 1
#define __FLT_MAX_10_EXP__ 38
#define __LONG_MAX__ 2147483647L
#define __DEC128_SUBNORMAL_MIN__ 0.000000000000000000000000000000001E-6143DL
#define __FLT_HAS_INFINITY__ 1
#define __DEC64_MAX__ 9.999999999999999E384DD
#define __CHAR16_TYPE__ short unsigned int
#define __DEC64_MANT_DIG__ 16
#define __DEC32_MAX_EXP__ 97
#define linux 1
#define __EXCEPTIONS 1
#define __LDBL_MANT_DIG__ 64
#define __DBL_HAS_QUIET_NAN__ 1
#define __WCHAR_TYPE__ int
#define __SIZEOF_FLOAT__ 4
#define __DEC64_MIN_EXP__ (-382)
#define __FLT_DIG__ 6
#define __INT_MAX__ 2147483647
#define __i486 1
#define __FLT_MAX_EXP__ 128
#define __DBL_MANT_DIG__ 53
#define __DEC64_MIN__ 1E-383DD
#define __WINT_TYPE__ unsigned int
#define __SIZEOF_SHORT__ 2
#define __LDBL_MIN_EXP__ (-16381)
#define __SSP__ 1
#define __LDBL_MAX_10_EXP__ 4932
#define __DBL_EPSILON__ 2.2204460492503131e-16
#define __SIZEOF_WCHAR_T__ 4
#define __DEC_EVAL_METHOD__ 2
#define __INTMAX_MAX__ 9223372036854775807LL
#define __FLT_DENORM_MIN__ 1.40129846e-45F
#define __CHAR32_TYPE__ unsigned int
#define __FLT_MAX__ 3.40282347e+38F
#define __SIZEOF_DOUBLE__ 8
#define __INTMAX_TYPE__ long long int
#define i386 1
#define __DEC128_MAX_EXP__ 6145
#define __DEC32_MANT_DIG__ 7
#define __DBL_MAX_10_EXP__ 308
#define __LDBL_DENORM_MIN__ 3.64519953188247460253e-4951L
#define __STDC__ 1
#define __PTRDIFF_TYPE__ int
#define __DEC64_SUBNORMAL_MIN__ 0.000000000000001E-383DD
#define __DEC128_MANT_DIG__ 34
#define __LDBL_MIN_10_EXP__ (-4931)
#define __SIZEOF_LONG_LONG__ 8
#define __LDBL_DIG__ 18
#define __GNUC_GNU_INLINE__ 1
#define _GNU_SOURCE 1

/* DQ (10/22/2010): Isolate some of the OSX specific details. */
#ifdef __APPLE__
/* DQ (10/13/2010): This is a specific fix suggested by Cong Hou for an OSX problem. 
 * The use of __BLOCKS__ (being defined as a macro) controls the code that ROSE will 
 * see in the Mac OSX header files.  If it is defined then ROSE will try to parse 
 * some nonstandard code that represents a Max OSX (or Apple) language extension.
 * So in this fix we detect if we are on an OSX OS (thus using OSX header files)
 * and we undefine the __BLOCKS__ macro.
 */
#undef __BLOCKS__

/* Cong (10/20/2010): In EDG 3.x, __attribute__, __dead2, __pure2 cannot be parsed, 
 * so we define them to nothing here.  */
#ifndef ROSE_USE_EDG_VERSION_4
	#define __attribute__(x)
	#define __dead2
	#define __pure2
#endif

#endif //#ifdef __APPLE__

/* This does not work since the include path may not be set to include this file. 
   Also including this is problematic since it would be used in both the compilation 
   of ROSE source code and source code using ROSE (which is awkward).
 */
/* Include rose_config.h (so that we can see and use the macro definitions) */
/* #include "rose_config.h" */

/* DQ (7/11/2009): fix emacs coloring... with a '$' */

/* Turn on use of restrict in EDG front-end using --edg:restrict */
#ifdef __GNUC__
/* for GNU g++ */
#define restrict __restrict__
#endif

/* Defined language modes in ROSE (permits control over what builting functions are defined) */
#define ROSE_C_LANGUAGE_MODE      0
#define ROSE_CXX_LANGUAGE_MODE    1
#define ROSE_CUDA_LANGUAGE_MODE   2
#define ROSE_OPENCL_LANGUAGE_MODE 3

#ifndef ROSE_LANGUAGE_MODE
  #error "Macro ROSE_LANGUAGE_MODE should have been defined ROSE_LANGUAGE_MODE == 0 for C and C99 and ROSE_LANGUAGE_MODE == 1 for C++ (CUDA == 2, OpenCL == 3)"
#endif

/* PC (9/25/2006): Define __builtin_va_start to the EDG expected symbol
   __builtin_stdarg_start */
#define __builtin_va_start __builtin_stdarg_start

/* Outside strict ISO C mode (-ansi, -std=c89 or -std=c99), the functions _exit, alloca,
    bcmp, bzero, dcgettext, dgettext, dremf, dreml, drem, exp10f, exp10l, exp10, ffsll,
    ffsl, ffs, fprintf_unlocked, fputs_unlocked, gammaf, gammal, gamma, gettext, index,
    isascii, j0f, j0l, j0, j1f, j1l, j1, jnf, jnl, jn, mempcpy, pow10f, pow10l, pow10,
    printf_unlocked, rindex, scalbf, scalbl, scalb, signbit, signbitf, signbitl,
    significandf, significandl, significand, sincosf, sincosl, sincos, stpcpy, strdup,
    strfmon, toascii, y0f, y0l, y0, y1f, y1l, y1, ynf, ynl and yn may be handled as
    built-in functions. All these functions have corresponding versions prefixed with
    __builtin_, which may be used even in strict C89 mode.
*/
/* extern void _exit ( int Status); */

/* The ISO C99 functions _Exit, acoshf, acoshl, acosh, asinhf, asinhl, asinh, atanhf,
    atanhl, atanh, cabsf, cabsl, cabs, cacosf, cacoshf, cacoshl, cacosh, cacosl, cacos,
    cargf, cargl, carg, casinf, casinhf, casinhl, casinh, casinl, casin, catanf, catanhf,
    catanhl, catanh, catanl, catan, cbrtf, cbrtl, cbrt, ccosf, ccoshf, ccoshl, ccosh,
    ccosl, ccos, cexpf, cexpl, cexp, cimagf, cimagl, cimag, conjf, conjl, conj, copysignf,
    copysignl, copysign, cpowf, cpowl, cpow, cprojf, cprojl, cproj, crealf, creall, creal,
    csinf, csinhf, csinhl, csinh, csinl, csin, csqrtf, csqrtl, csqrt, ctanf, ctanhf,
    ctanhl, ctanh, ctanl, ctan, erfcf, erfcl, erfc, erff, erfl, erf, exp2f, exp2l, exp2,
    expm1f, expm1l, expm1, fdimf, fdiml, fdim, fmaf, fmal, fmaxf, fmaxl, fmax, fma, fminf,
    fminl, fmin, hypotf, hypotl, hypot, ilogbf, ilogbl, ilogb, imaxabs, isblank, iswblank,
    lgammaf, lgammal, lgamma, llabs, llrintf, llrintl, llrint, llroundf, llroundl,
    llround, log1pf, log1pl, log1p, log2f, log2l, log2, logbf, logbl, logb, lrintf,
    lrintl, lrint, lroundf, lroundl, lround, nearbyintf, nearbyintl, nearbyint,
    nextafterf, nextafterl, nextafter, nexttowardf, nexttowardl, nexttoward, remainderf,
    remainderl, remainder, remquof, remquol, remquo, rintf, rintl, rint, roundf, roundl,
    round, scalblnf, scalblnl, scalbln, scalbnf, scalbnl, scalbn, snprintf, tgammaf,
    tgammal, tgamma, truncf, truncl, trunc, vfscanf, vscanf, vsnprintf and vsscanf are
    handled as built-in functions except in strict ISO C90 mode (-ansi or -std=c89). 
*/
/* extern void _Exit ( int Status); */
double      __builtin_copysign  (double __builtin__x, double __builtin__y);
float       __builtin_copysignf (float __builtin__x, float __builtin__y);
long double __builtin_copysignl (long double __builtin__x, long double __builtin__y);

/* There are also built-in versions of the ISO C99 functions acosf, acosl, asinf, asinl,
    atan2f, atan2l, atanf, atanl, ceilf, ceill, cosf, coshf, coshl, cosl, expf, expl,
    fabsf, fabsl, floorf, floorl, fmodf, fmodl, frexpf, frexpl, ldexpf, ldexpl, log10f,
    log10l, logf, logl, modfl, modf, powf, powl, sinf, sinhf, sinhl, sinl, sqrtf, sqrtl,
    tanf, tanhf, tanhl and tanl that are recognized in any mode since ISO C90 reserves
    these names for the purpose to which ISO C99 puts them. All these functions have
    corresponding versions prefixed with __builtin_. 
*/

/* DQ (4/12/2005): Required to compile the gnu version 3.4.3 cmath header file. Defined
//                 as functions instead of macros to avoid constant propagation issues.
*/
float       __builtin_acosf  (float __builtin__x);
long double __builtin_acosl  (long double __builtin__x);
float       __builtin_asinf  (float __builtin__x);
long double __builtin_asinl  (long double __builtin__x);
float       __builtin_atanf  (float __builtin__x);
long double __builtin_atanl  (long double __builtin__x);
float       __builtin_atan2f (float __builtin__x,float __builtin__y);
long double __builtin_atan2l (long double __builtin__x,long double __builtin__y);
float       __builtin_ceilf  (float __builtin__x);
long double __builtin_ceill  (long double __builtin__x);
float       __builtin_coshf  (float __builtin__x);
long double __builtin_coshl  (long double __builtin__x);
float       __builtin_floorf (float __builtin__x);
long double __builtin_floorl (long double __builtin__x);
float       __builtin_fmodf  (float __builtin__x,float __builtin__y);
long double __builtin_fmodl  (long double __builtin__x,long double __builtin__y);
float       __builtin_frexpf (float __builtin__x,int *__builtin__y);
long double __builtin_frexpl (long double __builtin__x,int *__builtin__y);
float       __builtin_ldexpf (float __builtin__x,float __builtin__y);
long double __builtin_ldexpl (long double __builtin__x,long double __builtin__y);
float       __builtin_log10f (float __builtin__x);
long double __builtin_log10l (long double __builtin__x);
float       __builtin_modff  (float __builtin__x,float *__builtin__y);
long double __builtin_modfl  (long double __builtin__x,long double *__builtin__y);
float       __builtin_powf   (float __builtin__x,float __builtin__y);
long double __builtin_powl   (long double __builtin__x,long double __builtin__y);
float       __builtin_sinhf  (float __builtin__x);
long double __builtin_sinhl  (long double __builtin__x);
float       __builtin_tanf   (float __builtin__x);
long double __builtin_tanl   (long double __builtin__x);
float       __builtin_tanhf  (float __builtin__x);
long double __builtin_tanhl  (long double __builtin__x);

/* DQ (5/15/2006): Suggested by Christian Biesinger (working with Markus Schordan) */
long double __builtin_powil  (long double __builtin__x, int __builtin__i);
double      __builtin_powi   (double __builtin__x, int __builtin__i);
float       __builtin_powif  (float __builtin__x, int __builtin__i);


/* The ISO C94 functions iswalnum, iswalpha, iswcntrl, iswdigit, iswgraph, iswlower,
    iswprint, iswpunct, iswspace, iswupper, iswxdigit, towlower and towupper are handled
    as built-in functions except in strict ISO C90 mode (-ansi or -std=c89). 
*/

/* The ISO C90 functions abort, abs, acos, asin, atan2, atan, calloc, ceil, cosh, cos,
    exit, exp, fabs, floor, fmod, fprintf, fputs, frexp, fscanf, isalnum, isalpha,
    iscntrl, isdigit, isgraph, islower, isprint, ispunct, isspace, isupper, isxdigit,
    tolower, toupper, labs, ldexp, log10, log, malloc, memcmp, memcpy, memset, modf, pow,
    printf, putchar, puts, scanf, sinh, sin, snprintf, sprintf, sqrt, sscanf, strcat,
    strchr, strcmp, strcpy, strcspn, strlen, strncat, strncmp, strncpy, strpbrk, strrchr,
    strspn, strstr, tanh, tan, vfprintf, vprintf and vsprintf are all recognized as
    built-in functions unless -fno-builtin is specified (or -fno-builtin-function is
    specified for an individual function). All of these functions have corresponding
    versions prefixed with __builtin_.
*/


/* DQ (7/29/2005): declarations for builtin functions used by GNU, but 
// already defined in EDG (it seems that we can provide declarations for 
// them explicitly).  These should be marked as compiler generated in 
// the AST.
*/
char *      __builtin_strchr (const char *__builtin__s, int __builtin__c);
char *      __builtin_strrchr(const char *__builtin__s, int __builtin__c);
char *      __builtin_strpbrk(const char *__builtin__s, const char *__builtin__accept);
char *      __builtin_strstr (const char *__builtin__haystack, const char *__builtin__needle);
float       __builtin_nansf  (const char *__builtin__x);
double      __builtin_nans   (const char *__builtin__x);
long double __builtin_nansl  (const char *__builtin__x);
double      __builtin_fabs   (double      __builtin__x);
float       __builtin_fabsf  (float       __builtin__x);
long double __builtin_fabsl  (long double __builtin__x);
float       __builtin_cosf   (float       __builtin__x);
long double __builtin_cosl   (long double __builtin__x);
float       __builtin_sinf   (float       __builtin__x);
long double __builtin_sinl   (long double __builtin__x);
float       __builtin_sqrtf  (float       __builtin__x);
long double __builtin_sqrtl  (long double __builtin__x);

/* DQ (5/4/2010): Reported problem under Ubuntu 9.10 with GCC 4.4.1 and Boost 1.42.
*/
int __builtin_fpclassify (int, int, int, int, int, ...);

/* GCC provides built-in versions of the ISO C99 floating point comparison macros that
    avoid raising exceptions for unordered operands. They have the same names as the
    standard macros ( isgreater, isgreaterequal, isless, islessequal, islessgreater, and
    isunordered) , with __builtin_ prefixed. We intend for a library implementor to be
    able to simply #define each standard macro to its built-in equivalent.

DQ (6/19/2007): These might required math.h to be included.
*/
# ifndef __builtin_isgreater
#  define __builtin_isgreater(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x > __y; }))
# endif

/* Return nonzero value if X is greater than or equal to Y.  */
# ifndef __builtin_isgreaterequal
#  define __builtin_isgreaterequal(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x >= __y; }))
# endif

/* Return nonzero value if X is less than Y.  */
# ifndef __builtin_isless
#  define __builtin_isless(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x < __y; }))
# endif

/* Return nonzero value if X is less than or equal to Y.  */
# ifndef __builtin_islessequal
#  define __builtin_islessequal(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && __x <= __y; }))
# endif

/* Return nonzero value if either X is less than Y or Y is less than X.  */
# ifndef __builtin_islessgreater
#  define __builtin_islessgreater(x, y) \
  (__extension__							      \
   ({ __typeof__(x) __x = (x); __typeof__(y) __y = (y);			      \
      !isunordered (__x, __y) && (__x < __y || __y < __x); }))
# endif

/* Return nonzero value if arguments are unordered.  */
# ifndef __builtin_isunordered
#  define __builtin_isunordered(u, v) \
  (__extension__							      \
   ({ __typeof__(u) __u = (u); __typeof__(v) __v = (v);			      \
      fpclassify (__u) == FP_NAN || fpclassify (__v) == FP_NAN; }))
# endif

/* int __builtin_isgreater(x,y) isgreater(x,y);
int __builtin_isgreaterequal(x,y) isgreaterequal(x,y);
int __builtin_isless(x,y) isless(x,y);
int __builtin_islessequal(x,y) islessequal(x,y);
int __builtin_islessgreater(x,y) islessgreater(x,y);
int __builtin_isunordered(x,y) isunordered(x,y);
*/

/* GNU also supports a few other types of builtin functions: */
void * __builtin_return_address (unsigned int level);
void * __builtin_frame_address (unsigned int level);

/* Additional builtin functions that take or return types as arguments
(described at http://gcc.gnu.org/onlinedocs/gcc-4.0.3/gcc/Other-Builtins.html).
There are more difficult to reproduce except as macros that define them away:

int __builtin_types_compatible_p (type1, type2);
type __builtin_choose_expr (const_exp, exp1, exp2);
int __builtin_constant_p (exp);

DQ (6/19/2007): The definitions below are not correct, but should be portable, 
they will only be a problem is the resulting code is unparsed directly from 
the AST.  The detection of the use of these macros in the near future will help 
make this safer.
*/

/* This is not a correct test, but it is a weak form of equivalence that is portable */
#define __builtin_types_compatible_p(T1,T2) (sizeof(T1)==sizeof(T2))
/* This is not correct, but it should be portable */
#define __builtin_choose_expr(exp,T1,T2) (T1)
/* This is not correct, but it should be portable, make it always return false for now */
#define __builtin_constant_p(exp) (0)


/* Additional builtin functions:
(also from http://gcc.gnu.org/onlinedocs/gcc-4.0.3/gcc/Other-Builtins.html):
 */

long        __builtin_expect (long __builtin__exp, long __builtin__c);
void        __builtin_prefetch (const void *__builtin__addr, ...);
double      __builtin_huge_val (void);
float       __builtin_huge_valf (void);
long double __builtin_huge_vall (void);
double      __builtin_inf (void);
float       __builtin_inff (void);
long double __builtin_infl (void);
double      __builtin_nan (const char *__builtin__str);
float       __builtin_nanf (const char *__builtin__str);
long double __builtin_nanl (const char *__builtin__str);
double      __builtin_nans (const char *__builtin__str);
float       __builtin_nansf (const char *__builtin__str);
long double __builtin_nansl (const char *__builtin__str);
/* DQ (6/19/2007): Commented out because in interferred with existing 
                   builtin function defined as int __builtin_ffs (int __builtin__x); in EDG.
int __builtin_ffs (unsigned int __builtin__x); */
int __builtin_clz (unsigned int __builtin__x);
int __builtin_ctz (unsigned int __builtin__x);
int __builtin_popcount (unsigned int __builtin__x);
int __builtin_parity (unsigned int __builtin__x);
int __builtin_ffsl (unsigned long __builtin__x);
int __builtin_clzl (unsigned long __builtin__x);
int __builtin_ctzl (unsigned long __builtin__x);
int __builtin_popcountl (unsigned long __builtin__x);
int __builtin_parityl (unsigned long __builtin__x);
int __builtin_ffsll (unsigned long long __builtin__x);
int __builtin_clzll (unsigned long long __builtin__x);
int __builtin_ctzll (unsigned long long __builtin__x);
int __builtin_popcountll (unsigned long long __builtin__x);
int __builtin_parityll (unsigned long long __builtin__x);
double      __builtin_powi (double __builtin__x, int __builtin__y);
float       __builtin_powif (float __builtin__x, int __builtin__y);
long double __builtin_powil (long double __builtin__x, int __builtin__y);


/* DQ (6/19/2007): I think these defines can be eliminated now in favor of the 
                   builtin function prototypes.
   DQ (8/20/2006): Let the builtin values be equal to the largest possible values.
   DQ (5/20/2006): These should be defined to be appropriate values or defined as 
   function prototypes as others are below.  These can't be functions because they 
   are used to initialize static const variables.

#define __builtin_huge_valf() __FLT_MAX__
#define __builtin_nanf(string) 0
#define __builtin_huge_val() __DBL_MAX__
#define __builtin_nan(string) 0
#define __builtin_huge_vall() __LDBL_MAX__
#define __builtin_nanl(string) 0
*/

/* DQ (8/25/2009): Added another builtin function required to compile ROSE with ROSE. */
#ifdef _GLIBCXX_ATOMIC_BUILTINS
typedef int _Atomic_word;
_Atomic_word  __sync_fetch_and_add(volatile _Atomic_word* __mem, int __val);
#endif

/* DQ (8/25/2009): Added another builtin function required to compile ROSE with ROSE. (required for boost) */
/* #if ROSE_CPP_MODE */
#if (ROSE_LANGUAGE_MODE == ROSE_CXX_LANGUAGE_MODE)
/* These use C++ references and so can't be used in C mode (should not be used with EDG 4.3) */
#ifndef ROSE_USE_NEW_EDG_INTERFACE
int __sync_lock_test_and_set( int & v, int n );
int __sync_lock_release( int & v );
#endif
#endif

/* DQ and Liao (7/11/2009) Added macros to define away new GNU C++ extension (not required for newer EDG 4.0 use */
// Only for GCC 4.3.0 and above
#if __GNUC__ > 4 || \
  (__GNUC__ == 4 && (__GNUC_MINOR__ > 3 || \
		     (__GNUC_MINOR__ == 3 && \
		      __GNUC_PATCHLEVEL__ >= 0)))

#define __has_nothrow_assign sizeof
#define __has_nothrow_copy sizeof
#define __has_nothrow_constructor sizeof
#define __has_trivial_assign sizeof
#define __has_trivial_copy  sizeof
#define __has_trivial_constructor sizeof
#define __has_trivial_destructor sizeof
#define __has_virtual_destructor sizeof
#define __is_abstract sizeof
//#define __is_base_of (base_type, derived_type)
#define __is_class sizeof
#define __is_empty sizeof
#define __is_enum sizeof
#define __is_pod sizeof 
#define __is_polymorphic sizeof
#define __is_union sizeof
void* __builtin_memmove(void * target, const void * source, unsigned long long nBytes);
void* __builtin_memchr(const  void * ptr, int value, unsigned long long num);

/* DQ (9/30/2009): This needs to be unsigned long long for 64-bit and 
   unsigned int for 32-bit (at least on GNU g++ versions 4.3 and greater).
   I think it is better to use size_t than __SIZE_TYPE__. 
   Note that this is a fix for a problem on 32-bit gnu 4.3 systems (nmi:x86_rhas_4).
// void* __builtin_memcpy (void * destination, const void * source, unsigned long long num );
// int __builtin_memcmp ( const void * ptr1, const void * ptr2, unsigned long long num );
// void* __builtin_memcpy (void * destination, const void * source, __SIZE_TYPE__ num );
// int __builtin_memcmp ( const void * ptr1, const void * ptr2, __SIZE_TYPE__ num );
void* __builtin_memcpy (void * destination, const void * source, __SIZE_TYPE__ num );
int __builtin_memcmp ( const void * ptr1, const void * ptr2, __SIZE_TYPE__ num );
*/

// changed it in edg 3.3/src/sys_predef.c instead since va_list is not declared here
//int __builtin_vsnprintf(char *str, unsigned long long size, const char *format, va_list ap);

#endif
/*
Target specific builtin functions are available at:
  http://gcc.gnu.org/onlinedocs/gcc-4.0.3/gcc/Target-Builtins.html
    * Alpha Built-in Functions
    * ARM Built-in Functions
    * Blackfin Built-in Functions
    * FR-V Built-in Functions
    * X86 Built-in Functions
    * MIPS Paired-Single Support
    * PowerPC AltiVec Built-in Functions
    * SPARC VIS Built-in Functions 
*/

/*
Required builtin function as supported by Intel (gnu builtin functions supported by Intel
C++ compiler).  I am not sure if we should include anything special specific to this, it
appears to be a subset of the more complete handling above.

__builtin_abs
__builtin_labs
__builtin_cos
__builtin_cosf
__builtin_fabs
__builtin_fabsf
__builtin_memcmp
__builtin_memcpy
__builtin_sin
__builtin_sinf
__builtin_sqrt
__builtin_sqrtf
__builtin_strcmp
__builtin_strlen
__builtin_strncmp
__builtin_abort
__builtin_prefetch
__builtin_constant_p
__builtin_printf
__builtin_fprintf
__builtin_fscanf
__builtin_scanf
__builtin_fputs
__builtin_memset
__builtin_strcat
__builtin_strcpy
__builtin_strncpy
__builtin_exit
__builtin_strchr
__builtin_strspn
__builtin_strcspn
__builtin_strstr
__builtin_strpbrk
__builtin_strrchr
__builtin_strncat
__builtin_alloca
__builtin_ffs
__builtin_index
__builtin_rindex
__builtin_bcmp
__builtin_bzero
__builtin_sinl
__builtin_cosl
__builtin_sqrtl
__builtin_fabsl
__builtin_frame_address (IA-32 only)
__builtin_return_address (IA-32 only)
*/

/* DQ (6/19/2007): For handling of offsetof macro we can't build a function 
prototype so EDG provides a mechanism to support this (from the basics.h 
header file). However, we can use the one defined in the GNU header files is 
we only handle the __offsetof__ macro (so define it way).  Then the builtin
function is just defined as being the offsetof macro.
*/
#define __offsetof__
/*
// DQ (2/9/2008): Don't define __builtin_offsetof in terms of offsetof since this
// causes a recursively defined marco on Fedora Core relase 4, though it works 
// fine on Red Hat Enterprise release 9. New definition taken from common
// implementation of offset, but modified to address specific case of where C++
// defines a operator&() which can be a problem for the more common definition.
// This definition is used to defin the __builtin_offsetof
// #define __builtin_offsetof(t,memb) offsetof(t,memb)
// #define __builtin_offsetof(t,memb) ((size_t)&(((t *)0)->memb))
*/
#ifndef __cplusplus
#define __builtin_offsetof(TYPE, MEMBER) ((size_t) &((TYPE *)0)->MEMBER)
#else
/* The cast to "char &" below avoids problems with user-defined "operator &", which can appear in a POD type. */
#define __builtin_offsetof(TYPE, MEMBER)					\
  (__offsetof__ (reinterpret_cast <size_t>			\
                 (&reinterpret_cast <const volatile char &>	\
                  (static_cast<TYPE *> (0)->MEMBER))))
#endif

/* DQ (10/30/2005): Added to allow compilation of g++ 
complex header, but this is likely the wrong thing to do.
This allows us to compile C++ code that uses the complex header file
but it turns "__complex__ double x" into just "double x". All this
is because the "__complex__" keyword is not supported in EDG.  Not
too much of a problem except that it means that the member function:
    complex(_ComplexT __z) : _M_value(__z) { }
needs to be specified using the "explicit" keyword.   Note that
none of this is required in g++ or icc, just in our version of 
ROSE using the g++ header files (because __complex__ is not
defined in EDG and the only thing to map it to is "", I think).

The following solution does NOT work.  Replacing
    typedef __complex__ double _ComplexT;
with
    typedef complex<double> _ComplexT;
in the complex header file.

The following solution does work.  Replacing
    complex(_ComplexT __z) : _M_value(__z) { }
with
    explicit complex(_ComplexT __z) : _M_value(__z) { }
in the complex header file.  It is not a great solution
and it would have to be done in the configuration of ROSE
to the complex header file.  It would also fool any analysis
of the complex class into thinking that the internal type 
was just float or double for complex<float> or complex<double>.
 */

/* DQ (8/22/2006):
   EDG does not appear to support __complex__ as a keyword.
   Thus we have to define it to white space to permit codes using 
   __complex__ to be compiled.  This appears to work fine, but
   we are likely to confuse an analysis that depends upon recognizing 
   complex types.

   Note that _Complex is a C99 type, also fequently recognised by C 
   compilers, so for non C++ codes we can translate __complex__ to 
   _Complex
*/
/* #ifdef __cplusplus */
/* #if ROSE_CPP_MODE */
#if (ROSE_LANGUAGE_MODE == ROSE_CXX_LANGUAGE_MODE)
/* C++ case is different because the header files will define a complex type (see hdrs1/complex) */

  #define __complex__
  #define __real__ 
  #define __imag__

/* DQ (9/26/2006): Plum Hall uses this for C++ code, but I think that we have define it away */
  #define _Complex
#else
/* This works for both C and C99 modes */
  #define __complex__ _Complex
  #define __real__ 
  #define __imag__
#endif

/* gcc uses the C99 name _Complex_I in <complex.h>, but our EDG doesn't handle
   the GCC extension that they define it to. */
#ifdef __INTEL_COMPILER
/* #define _Complex_I __I__ */
#else
#define _Complex_I __I__
#endif

/* Disable inclusion of complex.h on Linux */
#define _COMPLEX_H
/* Disable inclusion of complex.h on Mac OS X */
#define __COMPLEX__

/* DQ (8/25/2009): Added complex builtin functions as required to compile parts of ROSE with ROSE */
/* #if _GLIBCXX_USE_C99_COMPLEX */
/* #if _GLIBCXX_USE_C99_COMPLEX || defined(USE_ROSE) */
#if _GLIBCXX_USE_C99_COMPLEX
float  __builtin_cabsf(__complex__ float __z);
double __builtin_cabs (__complex__ double __z);
long double __builtin_cabsl(const __complex__ long double& __z);

float  __builtin_cargf(__complex__ float __z);
double __builtin_carg(__complex__ double __z);
long double __builtin_cargl(const __complex__ long double& __z);

__complex__ float  __builtin_ccosf(__complex__ float __z);
__complex__ double __builtin_ccos(__complex__ double __z);
__complex__ long double __builtin_ccosl(const __complex__ long double& __z);

__complex__ float  __builtin_ccoshf(__complex__ float __z);
__complex__ double __builtin_ccosh(__complex__ double __z);
__complex__ long double __builtin_ccoshl(const __complex__ long double& __z);

__complex__ float  __builtin_cexpf(__complex__ float __z);
__complex__ double __builtin_cexp(__complex__ double __z);
__complex__ long double __builtin_cexpl(const __complex__ long double& __z);

__complex__ float  __builtin_clogf(__complex__ float __z);
__complex__ double __builtin_clog(__complex__ double __z);
__complex__ long double __builtin_clogl(const __complex__ long double& __z);

__complex__ float  __builtin_csinf(__complex__ float __z);
__complex__ double __builtin_csin(__complex__ double __z);
__complex__ long double __builtin_csinl(const __complex__ long double& __z);

__complex__ float  __builtin_csinhf(__complex__ float __z);
__complex__ double __builtin_csinh(__complex__ double __z);
__complex__ long double __builtin_csinhl(const __complex__ long double& __z);

__complex__ float __builtin_csqrtf(__complex__ float __z);
__complex__ double __builtin_csqrt(__complex__ double __z);
__complex__ long double __builtin_csqrtl(const __complex__ long double& __z);

__complex__ float  __builtin_ctanf (__complex__ float __x);
__complex__ double __builtin_ctan(__complex__ double __z);
__complex__ long double __builtin_ctanl(const __complex__ long double& __z);

__complex__ float  __builtin_ctanhf(__complex__ float __z);
__complex__ double __builtin_ctanh(__complex__ double __z);
__complex__ long double __builtin_ctanhl(const __complex__ long double& __z);

__complex__ float  __builtin_cpowf(__complex__ float __x, __complex__ float __y);
__complex__ double __builtin_cpow(__complex__ double __x, __complex__ double __y);
__complex__ long double __builtin_cpowl(const __complex__ long double& __x,const __complex__ long double& __y);
#endif

#ifdef USE_ROSE
/* DQ (1/26/2010): Define these so that ROSE will compile with ROSE. */
#define __builtin_cabsf(x) 0
#define __builtin_cabs(x)  0
#define __builtin_cabsl(x) 0
#define __builtin_cargf(x) 0
#define __builtin_carg(x) 0
#define __builtin_cargl(x) 0
#define __builtin_ccosf(x) 0
#define __builtin_ccos(x) 0
#define __builtin_ccosl(x) 0
#define __builtin_ccoshf(x) 0
#define __builtin_ccosh(x) 0
#define __builtin_ccoshl(x) 0
#define __builtin_cexpf(x) 0
#define __builtin_cexp(x) 0
#define __builtin_cexpl(x) 0
#define __builtin_clogf(x) 0
#define __builtin_clog(x) 0
#define __builtin_clogl(x) 0
#define __builtin_csinf(x) 0
#define __builtin_csin(x) 0
#define __builtin_csinl(x) 0
#define __builtin_csinhf(x) 0
#define __builtin_csinh(x) 0
#define __builtin_csinhl(x) 0
#define __builtin_csqrtf(x) 0
#define __builtin_csqrt(x) 0
#define __builtin_csqrtl(x) 0
#define __builtin_ctanf(x) 0
#define __builtin_ctan(x) 0
#define __builtin_ctanl(x) 0
#define __builtin_ctanhf(x) 0
#define __builtin_ctanh(x) 0
#define __builtin_ctanhl(x) 0
#define __builtin_cpowf(x,y) 0
#define __builtin_cpow(x,y) 0
#define __builtin_cpowl(x,y) 0
#endif


/* Defined this to avoid warnings (e.g. test2001_11.C) from 3.4.6 systems header files. */
#define __weakref__(NAME)

/* DQ (6/19/2007): I think we can comment this out now, since it is better defined above!
   DQ (1/31/2007): GNU modifier required to handle code using the offsetof macro in C++ g++ 3.4 and greater */
/* #define __offsetof__ */

#if 0
/* DQ (7/11/2009) Added macros to define away new GNU C++ extension (not required for newer EDG 4.0 use */
#define __has_nothrow_assign (type)
#define __has_nothrow_copy (type)
#define __has_nothrow_constructor (type)
#define __has_trivial_assign (type)
#define __has_trivial_copy (type)
#define __has_trivial_constructor (type)
#define __has_trivial_destructor (type)
#define __has_virtual_destructor (type)
#define __is_abstract (type)
#define __is_base_of (base_type, derived_type)
#define __is_class (type)
#define __is_empty (type)
#define __is_enum (type)
#define __is_pod (type)
#define __is_polymorphic (type)
#define __is_union (type)
#endif

/* 
   DQ (7/15/2009): Added support for MS Windows Code 
   It might be that this file in included too late when using WINE.
   DQ (8/15/2009): The problem is that the rose_config.h file is 
   not included so USE_ROSE_WINDOWS_ANALYSIS_SUPPORT is not defined.
 */
/* #ifdef USE_ROSE_WINDOWS_ANALYSIS_SUPPORT */
 #define __builtin_ms_va_list __builtin_va_list
 #define __ms_va_list va_list
/* #endif */

/*
   DQ (9/12/2009): Avoid this GNU extension since it is a problem for EDG.
 */
#define _GLIBCXX_EXTERN_TEMPLATE 0

/*************************************************
        CUDA and OPENCL SPECIFIC SUPPORT
*************************************************/

#if (ROSE_LANGUAGE_MODE == ROSE_CUDA_LANGUAGE_MODE)

/* CUDA Built-in Types */

  /*
    Extract and adapted from:
      CUDA_PATH/include/common_functions.h
      CUDA_PATH/include/common_types.h
      CUDA_PATH/include/device_functions.h
      CUDA_PATH/include/device_types.h
      CUDA_PATH/include/sm_11_atomic_functions.h
      CUDA_PATH/include/sm_12_atomic_functions.h
      CUDA_PATH/include/sm_13_atomic_functions.h
      CUDA_PATH/include/sm_20_atomic_functions.h
      CUDA_PATH/include/sm_20_intrinsics.h
      CUDA_PATH/include/vector_types.h
      CUDA_PATH/include/vector_functions.h
  */

  /* Vector Types */

struct char1
{
  signed char x;
};

struct uchar1 
{
  unsigned char x;
};

struct __builtin_align__(2) char2
{
  signed char x, y;
};

struct __builtin_align__(2) uchar2
{
  unsigned char x, y;
};

struct char3
{
  signed char x, y, z;
};

struct uchar3
{
  unsigned char x, y, z;
};

struct __builtin_align__(4) char4
{
  signed char x, y, z, w;
};

struct __builtin_align__(4) uchar4
{
  unsigned char x, y, z, w;
};

struct short1
{
  short x;
};

struct ushort1
{
  unsigned short x;
};

struct __builtin_align__(4) short2
{
  short x, y;
};

struct __builtin_align__(4) ushort2
{
  unsigned short x, y;
};

struct short3
{
  short x, y, z;
};

struct ushort3
{
  unsigned short x, y, z;
};

struct __builtin_align__(8) short4
{
  short x, y, z, w;
};

struct __builtin_align__(8) ushort4
{
  unsigned short x, y, z, w;
};

struct int1
{
  int x;
};

struct uint1
{
  unsigned int x;
};

struct __builtin_align__(8) int2
{
  int x, y;
};

struct __builtin_align__(8) uint2
{
  unsigned int x, y;
};

struct int3
{
  int x, y, z;
};

struct uint3
{
  unsigned int x, y, z;
};

struct __builtin_align__(16) int4
{
  int x, y, z, w;
};

struct __builtin_align__(16) uint4
{
  unsigned int x, y, z, w;
};

struct long1
{
  long int x;
};

struct ulong1
{
  unsigned long x;
};

struct 
#if defined (_WIN32)
       __builtin_align__(8)
#else /* _WIN32 */
       __builtin_align__(2*sizeof(long int))
#endif /* _WIN32 */
                                             long2
{
  long int x, y;
};

struct 
#if defined (_WIN32)
       __builtin_align__(8)
#else /* _WIN32 */
       __builtin_align__(2*sizeof(unsigned long int))
#endif /* _WIN32 */
                                                      ulong2
{
  unsigned long int x, y;
};

#if !defined(__LP64__)

struct long3
{
  long int x, y, z;
};

struct ulong3
{
  unsigned long int x, y, z;
};

struct __builtin_align__(16) long4
{
  long int x, y, z, w;
};

struct __builtin_align__(16) ulong4
{
  unsigned long int x, y, z, w;
};

#endif /* !__LP64__ */

struct float1
{
  float x;
};

struct __builtin_align__(8) float2
{
  float x, y;
};

struct float3
{
  float x, y, z;
};

struct __builtin_align__(16) float4
{
  float x, y, z, w;
};

struct longlong1
{
  long long int x;
};

struct ulonglong1
{
  unsigned long long int x;
};

struct __builtin_align__(16) longlong2
{
  long long int x, y;
};

struct __builtin_align__(16) ulonglong2
{
  unsigned long long int x, y;
};

struct double1
{
  double x;
};

struct __builtin_align__(16) double2
{
  double x, y;
};

typedef struct char1 char1;
typedef struct uchar1 uchar1;
typedef struct char2 char2;
typedef struct uchar2 uchar2;
typedef struct char3 char3;
typedef struct uchar3 uchar3;
typedef struct char4 char4;
typedef struct uchar4 uchar4;
typedef struct short1 short1;
typedef struct ushort1 ushort1;
typedef struct short2 short2;
typedef struct ushort2 ushort2;
typedef struct short3 short3;
typedef struct ushort3 ushort3;
typedef struct short4 short4;
typedef struct ushort4 ushort4;
typedef struct int1 int1;
typedef struct uint1 uint1;
typedef struct int2 int2;
typedef struct uint2 uint2;
typedef struct int3 int3;
typedef struct uint3 uint3;
typedef struct int4 int4;
typedef struct uint4 uint4;
typedef struct long1 long1;
typedef struct ulong1 ulong1;
typedef struct long2 long2;
typedef struct ulong2 ulong2;
typedef struct long3 long3;
typedef struct ulong3 ulong3;
typedef struct long4 long4;
typedef struct ulong4 ulong4;
typedef struct float1 float1;
typedef struct float2 float2;
typedef struct float3 float3;
typedef struct float4 float4;
typedef struct longlong1 longlong1;
typedef struct ulonglong1 ulonglong1;
typedef struct longlong2 longlong2;
typedef struct ulonglong2 ulonglong2;
typedef struct double1 double1;
typedef struct double2 double2;

typedef struct dim3 dim3;

struct dim3
{
    unsigned int x, y, z;
#if defined(__cplusplus)
    dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1) : x(x), y(y), z(z) {}
    dim3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
    operator uint3(void) { uint3 t; t.x = x; t.y = y; t.z = z; return t; }
#endif /* __cplusplus */
};

/* CUDA Built-in Variables */

dim3  gridDim;
uint3 blockIdx;
dim3  blockDim;
uint3 threadIdx;
int   warpSize;


/* CUDA Built-in Functions */

  /* Vector Functions (constructors) */
  
static __inline__ __host__ __device__ char1 make_char1(signed char x);
static __inline__ __host__ __device__ uchar1 make_uchar1(unsigned char x);
static __inline__ __host__ __device__ char2 make_char2(signed char x, signed char y);
static __inline__ __host__ __device__ uchar2 make_uchar2(unsigned char x, unsigned char y);
static __inline__ __host__ __device__ char3 make_char3(signed char x, signed char y, signed char z);
static __inline__ __host__ __device__ uchar3 make_uchar3(unsigned char x, unsigned char y, unsigned char z);
static __inline__ __host__ __device__ char4 make_char4(signed char x, signed char y, signed char z, signed char w);
static __inline__ __host__ __device__ uchar4 make_uchar4(unsigned char x, unsigned char y, unsigned char z, unsigned char w);
static __inline__ __host__ __device__ short1 make_short1(short x);
static __inline__ __host__ __device__ ushort1 make_ushort1(unsigned short x);
static __inline__ __host__ __device__ short2 make_short2(short x, short y);
static __inline__ __host__ __device__ ushort2 make_ushort2(unsigned short x, unsigned short y);
static __inline__ __host__ __device__ short3 make_short3(short x,short y, short z);
static __inline__ __host__ __device__ ushort3 make_ushort3(unsigned short x, unsigned short y, unsigned short z);
static __inline__ __host__ __device__ short4 make_short4(short x, short y, short z, short w);
static __inline__ __host__ __device__ ushort4 make_ushort4(unsigned short x, unsigned short y, unsigned short z, unsigned short w);
static __inline__ __host__ __device__ int1 make_int1(int x);
static __inline__ __host__ __device__ uint1 make_uint1(unsigned int x);
static __inline__ __host__ __device__ int2 make_int2(int x, int y);
static __inline__ __host__ __device__ uint2 make_uint2(unsigned int x, unsigned int y);
static __inline__ __host__ __device__ int3 make_int3(int x, int y, int z);
static __inline__ __host__ __device__ uint3 make_uint3(unsigned int x, unsigned int y, unsigned int z);
static __inline__ __host__ __device__ int4 make_int4(int x, int y, int z, int w);
static __inline__ __host__ __device__ uint4 make_uint4(unsigned int x, unsigned int y, unsigned int z, unsigned int w);
static __inline__ __host__ __device__ long1 make_long1(long int x);
static __inline__ __host__ __device__ ulong1 make_ulong1(unsigned long int x);
static __inline__ __host__ __device__ long2 make_long2(long int x, long int y);
static __inline__ __host__ __device__ ulong2 make_ulong2(unsigned long int x, unsigned long int y);
static __inline__ __host__ __device__ long3 make_long3(long int x, long int y, long int z);
static __inline__ __host__ __device__ ulong3 make_ulong3(unsigned long int x, unsigned long int y, unsigned long int z);
static __inline__ __host__ __device__ long4 make_long4(long int x, long int y, long int z, long int w);
static __inline__ __host__ __device__ ulong4 make_ulong4(unsigned long int x, unsigned long int y, unsigned long int z, unsigned long int w);
static __inline__ __host__ __device__ float1 make_float1(float x);
static __inline__ __host__ __device__ float2 make_float2(float x, float y);
static __inline__ __host__ __device__ float3 make_float3(float x, float y, float z);
static __inline__ __host__ __device__ float4 make_float4(float x, float y, float z, float w);
static __inline__ __host__ __device__ longlong1 make_longlong1(long long int x);
static __inline__ __host__ __device__ ulonglong1 make_ulonglong1(unsigned long long int x);
static __inline__ __host__ __device__ longlong2 make_longlong2(long long int x, long long int y);
static __inline__ __host__ __device__ ulonglong2 make_ulonglong2(unsigned long long int x, unsigned long long int y);
static __inline__ __host__ __device__ double1 make_double1(double x);
static __inline__ __host__ __device__ double2 make_double2(double x, double y);
  
  /* Synchronization functions */

__device__ void __threadfence_block();
__device__ void __threadfence();
__device__ void __threadfence_system();
__device__ void __syncthreads();
__device__ int  __syncthreads_count(int predicate);
__device__ int  __syncthreads_and(int predicate);
__device__ int  __syncthreads_or(int predicate);
  
  /* Time function */
  
__device__ clock_t clock();
  
  /* Atomic functions */

static __inline__ __device__ int atomicAdd(int *address, int val);
static __inline__ __device__ unsigned int atomicAdd(unsigned int *address, unsigned int val);
static __inline__ __device__ unsigned long long int atomicAdd(unsigned long long int *address, unsigned long long int val);
static __inline__ __device__ float atomicAdd(float *address, float val);
static __inline__ __device__ int atomicSub(int *address, int val);
static __inline__ __device__ unsigned int atomicSub(unsigned int *address, unsigned int val);
static __inline__ __device__ int atomicExch(int *address, int val);
static __inline__ __device__ unsigned int atomicExch(unsigned int *address, unsigned int val);
static __inline__ __device__ unsigned long long int atomicExch(unsigned long long int *address, unsigned long long int val);
static __inline__ __device__ float atomicExch(float *address, float val);
static __inline__ __device__ int atomicMin(int *address, int val);
static __inline__ __device__ unsigned int atomicMin(unsigned int *address, unsigned int val);
static __inline__ __device__ int atomicMax(int *address, int val);
static __inline__ __device__ unsigned int atomicMax(unsigned int *address, unsigned int val);
static __inline__ __device__ unsigned int atomicInc(unsigned int *address, unsigned int val);
static __inline__ __device__ unsigned int atomicDec(unsigned int *address, unsigned int val);
static __inline__ __device__ int atomicCAS(int *address, int compare, int val);
static __inline__ __device__ unsigned int atomicCAS(unsigned int *address, unsigned int compare, unsigned int val);
static __inline__ __device__ unsigned long long int atomicCAS(unsigned long long int *address, unsigned long long int compare, unsigned long long int val);
static __inline__ __device__ int atomicAnd(int *address, int val);
static __inline__ __device__ unsigned int atomicAnd(unsigned int *address, unsigned int val);
static __inline__ __device__ int atomicOr(int *address, int val);
static __inline__ __device__ unsigned int atomicOr(unsigned int *address, unsigned int val);
static __inline__ __device__ int atomicXor(int *address, int val);
static __inline__ __device__ unsigned int atomicXor(unsigned int *address, unsigned int val);

  /* Warp Vote functions */
  
__device__ int __all(int cond);
__device__ int __any(int cond);
__device__ unsigned int __ballot(int);

  /* Profiler Counter functions */

__device__ void __prof_trigger(int);

  /* Mathematical functions (TODO) */

  /* Texture functions (TODO) */

/* CUDA API (TODO) */

#endif


#if (ROSE_LANGUAGE_MODE == ROSE_OPENCL_LANGUAGE_MODE)

/* OpenCL Built-in Types */
  
  /* Extract and adapted from: CUDA_HOME/include/CL/cl_platform.h */

  /* Define alignment keys */

#if defined( __GNUC__ )
    #define OPENCL_ALIGNED(_x)          __attribute__ ((aligned(_x)))
#elif defined( _WIN32) && (_MSC_VER)
    /* Alignment keys neutered on windows because MSVC can't swallow function arguments with alignment requirements     */
    /* http://msdn.microsoft.com/en-us/library/373ak2y1%28VS.71%29.aspx                                                 */
    /* #include <crtdefs.h>                                                                                             */
    /* #define OPENCL_ALIGNED(_x)          _CRT_ALIGN(_x)                                                                   */
    #define OPENCL_ALIGNED(_x)
#else
   #warning  Need to implement some method to align data here
   #define  OPENCL_ALIGNED(_x)
#endif

  /* Scalar Types  */

#if (defined (_WIN32) && defined(_MSC_VER))

typedef signed   __int8         opencl_char;
typedef unsigned __int8         opencl_uchar;
typedef signed   __int16        opencl_short;
typedef unsigned __int16        opencl_ushort;
typedef signed   __int32        opencl_int;
typedef unsigned __int32        opencl_uint;
typedef signed   __int64        opencl_long;
typedef unsigned __int64        opencl_ulong;

typedef float                   opencl_float;
typedef double                  opencl_double;

typedef unsigned __int16        half;

typedef unsigned __int32        size_t;

#else

#include <stdint.h>

typedef int8_t          opencl_char;
typedef uint8_t         opencl_uchar;
typedef int16_t         opencl_short    __attribute__((aligned(2)));
typedef uint16_t        opencl_ushort   __attribute__((aligned(2)));
typedef int32_t         opencl_int      __attribute__((aligned(4)));
typedef uint32_t        opencl_uint     __attribute__((aligned(4)));
typedef int64_t         opencl_long     __attribute__((aligned(8)));
typedef uint64_t        opencl_ulong    __attribute__((aligned(8)));

typedef float           opencl_float    __attribute__((aligned(4)));
typedef double          opencl_double   __attribute__((aligned(8)));

typedef uint16_t        half     __attribute__((aligned(2)));

#endif

typedef uint32_t ptrdiff_t;
typedef uint32_t intptr_t;
typedef uint32_t uintptr_t;

  /* Vector Types */
  
typedef union
{
    opencl_char  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_char  x, y; };
   __extension__ struct { opencl_char  s0, s1; };
   __extension__ struct { opencl_char  lo, hi; };
#endif
} char2;

typedef union
{
    opencl_char  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_char  x, y, z, w; };
   __extension__ struct { opencl_char  s0, s1, s2, s3; };
   __extension__ struct { char2 lo, hi; };
#endif
} char4;

typedef union
{
    opencl_char   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_char  x, y, z, w; };
   __extension__ struct { opencl_char  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { char4 lo, hi; };
#endif
} char8;

typedef union
{
    opencl_char  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_char  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_char  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { char8 lo, hi; };
#endif
} char16;

typedef union
{
    opencl_uchar  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uchar  x, y; };
   __extension__ struct { opencl_uchar  s0, s1; };
   __extension__ struct { opencl_uchar  lo, hi; };
#endif
} uchar2;

typedef union
{
    opencl_uchar  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uchar  x, y, z, w; };
   __extension__ struct { opencl_uchar  s0, s1, s2, s3; };
   __extension__ struct { uchar2 lo, hi; };
#endif
} uchar4;

typedef union
{
    opencl_uchar   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uchar  x, y, z, w; };
   __extension__ struct { opencl_uchar  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { uchar4 lo, hi; };
#endif
} uchar8;

typedef union
{
    opencl_uchar  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uchar  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_uchar  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { uchar8 lo, hi; };
#endif
} uchar16;

typedef union
{
    opencl_short  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_short  x, y; };
   __extension__ struct { opencl_short  s0, s1; };
   __extension__ struct { opencl_short  lo, hi; };
#endif
} short2;

typedef union
{
    opencl_short  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_short  x, y, z, w; };
   __extension__ struct { opencl_short  s0, s1, s2, s3; };
   __extension__ struct { short2 lo, hi; };
#endif
} short4;

typedef union
{
    opencl_short   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_short  x, y, z, w; };
   __extension__ struct { opencl_short  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { short4 lo, hi; };
#endif
} short8;

typedef union
{
    opencl_short  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_short  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_short  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { short8 lo, hi; };
#endif
} short16;

typedef union
{
    opencl_ushort  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ushort  x, y; };
   __extension__ struct { opencl_ushort  s0, s1; };
   __extension__ struct { opencl_ushort  lo, hi; };
#endif
} ushort2;

typedef union
{
    opencl_ushort  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ushort  x, y, z, w; };
   __extension__ struct { opencl_ushort  s0, s1, s2, s3; };
   __extension__ struct { ushort2 lo, hi; };
#endif
} ushort4;

typedef union
{
    opencl_ushort   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ushort  x, y, z, w; };
   __extension__ struct { opencl_ushort  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { ushort4 lo, hi; };
#endif
} ushort8;

typedef union
{
    opencl_ushort  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ushort  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_ushort  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { ushort8 lo, hi; };
#endif
} ushort16;

typedef union
{
    opencl_int  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_int  x, y; };
   __extension__ struct { opencl_int  s0, s1; };
   __extension__ struct { opencl_int  lo, hi; };
#endif
} int2;

typedef union
{
    opencl_int  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_int  x, y, z, w; };
   __extension__ struct { opencl_int  s0, s1, s2, s3; };
   __extension__ struct { int2 lo, hi; };
#endif
} int4;

typedef union
{
    opencl_int   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_int  x, y, z, w; };
   __extension__ struct { opencl_int  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { int4 lo, hi; };
#endif
} int8;

typedef union
{
    opencl_int  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_int  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_int  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { int8 lo, hi; };
#endif
} int16;

typedef union
{
    opencl_uint  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uint  x, y; };
   __extension__ struct { opencl_uint  s0, s1; };
   __extension__ struct { opencl_uint  lo, hi; };
#endif
} uint2;

typedef union
{
    opencl_uint  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uint  x, y, z, w; };
   __extension__ struct { opencl_uint  s0, s1, s2, s3; };
   __extension__ struct { uint2 lo, hi; };
#endif
} uint4;

typedef union
{
    opencl_uint   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uint  x, y, z, w; };
   __extension__ struct { opencl_uint  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { uint4 lo, hi; };
#endif
} uint8;

typedef union
{
    opencl_uint  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_uint  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_uint  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { uint8 lo, hi; };
#endif
} uint16;

typedef union
{
    opencl_long  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_long  x, y; };
   __extension__ struct { opencl_long  s0, s1; };
   __extension__ struct { opencl_long  lo, hi; };
#endif
} long2;

typedef union
{
    opencl_long  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_long  x, y, z, w; };
   __extension__ struct { opencl_long  s0, s1, s2, s3; };
   __extension__ struct { long2 lo, hi; };
#endif
} long4;

typedef union
{
    opencl_long   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_long  x, y, z, w; };
   __extension__ struct { opencl_long  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { long4 lo, hi; };
#endif
} long8;

typedef union
{
    opencl_long  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_long  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_long  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { long8 lo, hi; };
#endif
} long16;

typedef union
{
    opencl_ulong  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ulong  x, y; };
   __extension__ struct { opencl_ulong  s0, s1; };
   __extension__ struct { opencl_ulong  lo, hi; };
#endif
} ulong2;

typedef union
{
    opencl_ulong  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ulong  x, y, z, w; };
   __extension__ struct { opencl_ulong  s0, s1, s2, s3; };
   __extension__ struct { ulong2 lo, hi; };
#endif
} ulong4;

typedef union
{
    opencl_ulong   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ulong  x, y, z, w; };
   __extension__ struct { opencl_ulong  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { ulong4 lo, hi; };
#endif
} ulong8;

typedef union
{
    opencl_ulong  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_ulong  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_ulong  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { ulong8 lo, hi; };
#endif
} ulong16;

typedef union
{
    opencl_float  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_float  x, y; };
   __extension__ struct { opencl_float  s0, s1; };
   __extension__ struct { opencl_float  lo, hi; };
#endif
} float2;

typedef union
{
    opencl_float  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_float  x, y, z, w; };
   __extension__ struct { opencl_float  s0, s1, s2, s3; };
   __extension__ struct { float2 lo, hi; };
#endif
} float4;

typedef union
{
    opencl_float   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_float  x, y, z, w; };
   __extension__ struct { opencl_float  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { float4 lo, hi; };
#endif
} float8;

typedef union
{
    opencl_float  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_float  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_float  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { float8 lo, hi; };
#endif
} float16;

typedef union
{
    opencl_double  OPENCL_ALIGNED(2) s[2];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_double  x, y; };
   __extension__ struct { opencl_double  s0, s1; };
   __extension__ struct { opencl_double  lo, hi; };
#endif
} double2;

typedef union
{
    opencl_double  OPENCL_ALIGNED(4) s[4];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_double  x, y, z, w; };
   __extension__ struct { opencl_double  s0, s1, s2, s3; };
   __extension__ struct { double2 lo, hi; };
#endif
} double4;

typedef union
{
    opencl_double   OPENCL_ALIGNED(8) s[8];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_double  x, y, z, w; };
   __extension__ struct { opencl_double  s0, s1, s2, s3, s4, s5, s6, s7; };
   __extension__ struct { double4 lo, hi; };
#endif
} double8;

typedef union
{
    opencl_double  OPENCL_ALIGNED(16) s[16];
#if defined( __GNUC__) && ! defined( __STRICT_ANSI__ )
   __extension__ struct { opencl_double  x, y, z, w, __spacer4, __spacer5, __spacer6, __spacer7, __spacer8, __spacer9, sa, sb, sc, sd, se, sf; };
   __extension__ struct { opencl_double  s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF; };
   __extension__ struct { double8 lo, hi; };
#endif
} double16;


/* OpenCL Built-in Functions */

  /* Work-Item Functions */

uint get_work_dim();
size_t get_global_size(uint dimindx);
size_t get_global_id(uint dimindx);
size_t get_local_size(uint dimindx);
size_t get_local_id(uint dimindx);
size_t get_num_groups(uint dimindx);
size_t get_group_id(uint dimindx);

  /* Image Functions */
  
float4 read_imagef  (image2d_t image, sampler_t sampler, int2 coord);
float4 read_imagef  (image2d_t image, sampler_t sampler, float2 coord);
int4   read_imagei  (image2d_t image, sampler_t sampler, int2 coord);
int4   read_imagei  (image2d_t image, sampler_t sampler, float2 coord);
uint4  read_imageui (image2d_t image, sampler_t sampler, int2 coord);
uint4  read_imageui (image2d_t image, sampler_t sampler, float2 coord);
void   write_imagef (image2d_t image, int2 coord, ﬂoat4 color);
void   write_imagei (image2d_t image, int2 coord, int4 color);
void   write_imageui(image2d_t image, int2 coord, uint4 color);
float4 read_imagef  (image3d_t image, sampler_t sampler, int4 coord);
float4 read_imagef  (image3d_t image, sampler_t sampler, float4 coord);
int4   read_imagei  (image3d_t image, sampler_t sampler, int4 coord);
int4   read_imagei  (image3d_t image, sampler_t sampler, float4 coord);
uint4  read_imageui (image3d_t image, sampler_t sampler, int4 coord);
uint4  read_imageui (image3d_t image, sampler_t sampler, float4 coord);
int    get_image_width (image2d_t image);
int    get_image_width (image3d_t image);
int    get_image_height(image2d_t image);
int    get_image_height(image3d_t image);
int    get_image_depht (image3d_t image);
int    get_image_channel_data_type (image2d_t image);
int    get_image_channel_data_type (image3d_t image);
int    get_image_channel_order (image2d_t image);
int    get_image_channel_order (image3d_t image);
int2   get_image_dim (image2d_t image);
int4   get_image_dim (image3d_t image);

  /* Synchronization Functions */
  
void barrier (cl_mem_fence_ﬂags ﬂags);

  /* Explicit Memory Fence Functions */
  
void mem_fence (cl_mem_fence_ﬂags ﬂags);
void read_mem_fence (cl_mem_fence_ﬂags ﬂags);
void write_mem_fence (cl_mem_fence_ﬂags ﬂags);

  /* Miscellaneous Functions */
  
void wait_group_events (int num_events, event_t * event_list);
event_t async_work_group_copy (__local char *dst, const __global char *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global char *dst, const __local char *src, size_t num_elements, event_t event);
void prefetch (const __global char *p, size_t num_elements);
event_t async_work_group_copy (__local char2 *dst, const __global char2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global char2 *dst, const __local char2 *src, size_t num_elements, event_t event);
void prefetch (const __global char2 *p, size_t num_elements);
event_t async_work_group_copy (__local char4 *dst, const __global char4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global char4 *dst, const __local char4 *src, size_t num_elements, event_t event);
void prefetch (const __global char4 *p, size_t num_elements);
event_t async_work_group_copy (__local char8 *dst, const __global char8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global char8 *dst, const __local char8 *src, size_t num_elements, event_t event);
void prefetch (const __global char8 *p, size_t num_elements);
event_t async_work_group_copy (__local char16 *dst, const __global char16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global char16 *dst, const __local char16 *src, size_t num_elements, event_t event);
void prefetch (const __global char16 *p, size_t num_elements);
event_t async_work_group_copy (__local uchar *dst, const __global uchar *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uchar *dst, const __local uchar *src, size_t num_elements, event_t event);
void prefetch (const __global uchar *p, size_t num_elements);
event_t async_work_group_copy (__local uchar2 *dst, const __global uchar2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uchar2 *dst, const __local uchar2 *src, size_t num_elements, event_t event);
void prefetch (const __global uchar2 *p, size_t num_elements);
event_t async_work_group_copy (__local uchar4 *dst, const __global uchar4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uchar4 *dst, const __local uchar4 *src, size_t num_elements, event_t event);
void prefetch (const __global uchar4 *p, size_t num_elements);
event_t async_work_group_copy (__local uchar8 *dst, const __global uchar8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uchar8 *dst, const __local uchar8 *src, size_t num_elements, event_t event);
void prefetch (const __global uchar8 *p, size_t num_elements);
event_t async_work_group_copy (__local uchar16 *dst, const __global uchar16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uchar16 *dst, const __local uchar16 *src, size_t num_elements, event_t event);
void prefetch (const __global uchar16 *p, size_t num_elements);
event_t async_work_group_copy (__local short *dst, const __global short *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global short *dst, const __local short *src, size_t num_elements, event_t event);
void prefetch (const __global short *p, size_t num_elements);
event_t async_work_group_copy (__local short2 *dst, const __global short2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global short2 *dst, const __local short2 *src, size_t num_elements, event_t event);
void prefetch (const __global short2 *p, size_t num_elements);
event_t async_work_group_copy (__local short4 *dst, const __global short4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global short4 *dst, const __local short4 *src, size_t num_elements, event_t event);
void prefetch (const __global short4 *p, size_t num_elements);
event_t async_work_group_copy (__local short8 *dst, const __global short8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global short8 *dst, const __local short8 *src, size_t num_elements, event_t event);
void prefetch (const __global short8 *p, size_t num_elements);
event_t async_work_group_copy (__local short16 *dst, const __global short16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global short16 *dst, const __local short16 *src, size_t num_elements, event_t event);
void prefetch (const __global short16 *p, size_t num_elements);
event_t async_work_group_copy (__local ushort *dst, const __global ushort *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ushort *dst, const __local ushort *src, size_t num_elements, event_t event);
void prefetch (const __global ushort *p, size_t num_elements);
event_t async_work_group_copy (__local ushort2 *dst, const __global ushort2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ushort2 *dst, const __local ushort2 *src, size_t num_elements, event_t event);
void prefetch (const __global ushort2 *p, size_t num_elements);
event_t async_work_group_copy (__local ushort4 *dst, const __global ushort4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ushort4 *dst, const __local ushort4 *src, size_t num_elements, event_t event);
void prefetch (const __global ushort4 *p, size_t num_elements);
event_t async_work_group_copy (__local ushort8 *dst, const __global ushort8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ushort8 *dst, const __local ushort8 *src, size_t num_elements, event_t event);
void prefetch (const __global ushort8 *p, size_t num_elements);
event_t async_work_group_copy (__local ushort16 *dst, const __global ushort16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ushort16 *dst, const __local ushort16 *src, size_t num_elements, event_t event);
void prefetch (const __global ushort16 *p, size_t num_elements);
event_t async_work_group_copy (__local int *dst, const __global int *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global int *dst, const __local int *src, size_t num_elements, event_t event);
void prefetch (const __global int *p, size_t num_elements);
event_t async_work_group_copy (__local int2 *dst, const __global int2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global int2 *dst, const __local int2 *src, size_t num_elements, event_t event);
void prefetch (const __global int2 *p, size_t num_elements);
event_t async_work_group_copy (__local int4 *dst, const __global int4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global int4 *dst, const __local int4 *src, size_t num_elements, event_t event);
void prefetch (const __global int4 *p, size_t num_elements);
event_t async_work_group_copy (__local int8 *dst, const __global int8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global int8 *dst, const __local int8 *src, size_t num_elements, event_t event);
void prefetch (const __global int8 *p, size_t num_elements);
event_t async_work_group_copy (__local int16 *dst, const __global int16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global int16 *dst, const __local int16 *src, size_t num_elements, event_t event);
void prefetch (const __global int16 *p, size_t num_elements);
event_t async_work_group_copy (__local uint *dst, const __global uint *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uint *dst, const __local uint *src, size_t num_elements, event_t event);
void prefetch (const __global uint *p, size_t num_elements);
event_t async_work_group_copy (__local uint2 *dst, const __global uint2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uint2 *dst, const __local uint2 *src, size_t num_elements, event_t event);
void prefetch (const __global uint2 *p, size_t num_elements);
event_t async_work_group_copy (__local uint4 *dst, const __global uint4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uint4 *dst, const __local uint4 *src, size_t num_elements, event_t event);
void prefetch (const __global uint4 *p, size_t num_elements);
event_t async_work_group_copy (__local uint8 *dst, const __global uint8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uint8 *dst, const __local uint8 *src, size_t num_elements, event_t event);
void prefetch (const __global uint8 *p, size_t num_elements);
event_t async_work_group_copy (__local uint16 *dst, const __global uint16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global uint16 *dst, const __local uint16 *src, size_t num_elements, event_t event);
void prefetch (const __global uint16 *p, size_t num_elements);
event_t async_work_group_copy (__local long *dst, const __global long *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global long *dst, const __local long *src, size_t num_elements, event_t event);
void prefetch (const __global long *p, size_t num_elements);
event_t async_work_group_copy (__local long2 *dst, const __global long2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global long2 *dst, const __local long2 *src, size_t num_elements, event_t event);
void prefetch (const __global long2 *p, size_t num_elements);
event_t async_work_group_copy (__local long4 *dst, const __global long4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global long4 *dst, const __local long4 *src, size_t num_elements, event_t event);
void prefetch (const __global long4 *p, size_t num_elements);
event_t async_work_group_copy (__local long8 *dst, const __global long8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global long8 *dst, const __local long8 *src, size_t num_elements, event_t event);
void prefetch (const __global long8 *p, size_t num_elements);
event_t async_work_group_copy (__local long16 *dst, const __global long16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global long16 *dst, const __local long16 *src, size_t num_elements, event_t event);
void prefetch (const __global long16 *p, size_t num_elements);
event_t async_work_group_copy (__local ulong *dst, const __global ulong *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ulong *dst, const __local ulong *src, size_t num_elements, event_t event);
void prefetch (const __global ulong *p, size_t num_elements);
event_t async_work_group_copy (__local ulong2 *dst, const __global ulong2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ulong2 *dst, const __local ulong2 *src, size_t num_elements, event_t event);
void prefetch (const __global ulong2 *p, size_t num_elements);
event_t async_work_group_copy (__local ulong4 *dst, const __global ulong4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ulong4 *dst, const __local ulong4 *src, size_t num_elements, event_t event);
void prefetch (const __global ulong4 *p, size_t num_elements);
event_t async_work_group_copy (__local ulong8 *dst, const __global ulong8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ulong8 *dst, const __local ulong8 *src, size_t num_elements, event_t event);
void prefetch (const __global ulong8 *p, size_t num_elements);
event_t async_work_group_copy (__local ulong16 *dst, const __global ulong16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global ulong16 *dst, const __local ulong16 *src, size_t num_elements, event_t event);
void prefetch (const __global ulong16 *p, size_t num_elements);
event_t async_work_group_copy (__local float *dst, const __global float *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global float *dst, const __local float *src, size_t num_elements, event_t event);
void prefetch (const __global float *p, size_t num_elements);
event_t async_work_group_copy (__local float2 *dst, const __global float2 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global float2 *dst, const __local float2 *src, size_t num_elements, event_t event);
void prefetch (const __global float2 *p, size_t num_elements);
event_t async_work_group_copy (__local float4 *dst, const __global float4 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global float4 *dst, const __local float4 *src, size_t num_elements, event_t event);
void prefetch (const __global float4 *p, size_t num_elements);
event_t async_work_group_copy (__local float8 *dst, const __global float8 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global float8 *dst, const __local float8 *src, size_t num_elements, event_t event);
void prefetch (const __global float8 *p, size_t num_elements);
event_t async_work_group_copy (__local float16 *dst, const __global float16 *src, size_t num_elements, event_t event);
event_t async_work_group_copy (__global float16 *dst, const __local float16 *src, size_t num_elements, event_t event);
void prefetch (const __global float16 *p, size_t num_elements);

  /* Miscellaneous Functions (TODO) */



#endif


/***************************************************************************************
  SSE and MMX support (declarations for builtin functions defined in GNU header files)
****************************************************************************************/

/* From: mmintrin.h */
void __builtin_ia32_emms();
int __builtin_ia32_vec_init_v2si(int,int);
int __builtin_ia32_vec_ext_v2si(int,int);
int __builtin_ia32_packsswb(short,short);
int __builtin_ia32_packssdw(int,int);
int __builtin_ia32_packuswb(short,short);
int __builtin_ia32_punpckhbw(short,short);
int __builtin_ia32_punpckhwd(short,short);
int __builtin_ia32_punpckhdq(int,int);
int __builtin_ia32_punpcklbw(char,char);
int __builtin_ia32_punpcklwd(short,short);
int __builtin_ia32_punpckldq(int,int);
int __builtin_ia32_paddb(char,char);
int __builtin_ia32_paddw(int,int);
int __builtin_ia32_paddd(short,short);
int __builtin_ia32_paddq(long long,long long);
int __builtin_ia32_paddsb(char,char);
int __builtin_ia32_paddsw(int,int);
int __builtin_ia32_paddusb(char,char);
int __builtin_ia32_paddusw(int,int);
int __builtin_ia32_psubb(char,char);
int __builtin_ia32_psubw(int,int);
int __builtin_ia32_psubd(int,int);
int __builtin_ia32_psubq(long long,long long);
int __builtin_ia32_psubsb(char,char);
int __builtin_ia32_psubsw(short,short);
int __builtin_ia32_psubusb(char,char);
int __builtin_ia32_psubusw(int,int);
int __builtin_ia32_pmaddwd(short,short);
int __builtin_ia32_pmulhw(short,short);
int __builtin_ia32_pmullw(short,short);
int __builtin_ia32_psllw(short,long long);
int __builtin_ia32_pslld(int,long long);
int __builtin_ia32_psllq(long long, long long);
int __builtin_ia32_psraw(short, long long);
int __builtin_ia32_psrad(short, long long);
int __builtin_ia32_psrlw(short,long long);
int __builtin_ia32_psrld(short, long long);
int __builtin_ia32_psrlq(long long, long long);
int __builtin_ia32_pand(int,int);
int __builtin_ia32_pandn(int,int);
int __builtin_ia32_por(int,int);
int __builtin_ia32_pxor(int,int);
int __builtin_ia32_pcmpeqb(char,char);
int __builtin_ia32_pcmpgtb(char,char);
int __builtin_ia32_pcmpeqw(short,short);
int __builtin_ia32_pcmpgtw(short,short);
int __builtin_ia32_pcmpeqd(int,int);
int __builtin_ia32_pcmpgtd(int,int);
int __builtin_ia32_vec_init_v2si(int,int);
int __builtin_ia32_vec_init_v4hi(short,short,short,short);
int __builtin_ia32_vec_init_v8qi(char,char,char,char,char,char,char,char);

/* From: xmmintrin.h */
int __builtin_ia32_addss(float,float);
int __builtin_ia32_subss(float,float);
int __builtin_ia32_mulss(float,float);
int __builtin_ia32_divss(float,float);
int __builtin_ia32_sqrtss(float);
int __builtin_ia32_rcpss(float);
int __builtin_ia32_rsqrtss(float);
int __builtin_ia32_minss(float,float);
int __builtin_ia32_maxss(float,float);
int __builtin_ia32_addps(float,float);
int __builtin_ia32_subps(float,float);
int __builtin_ia32_mulps(float,float);
int __builtin_ia32_divps(float,float);
int __builtin_ia32_sqrtps(float);
int __builtin_ia32_rcpps(float);
int __builtin_ia32_rsqrtps(float);
int __builtin_ia32_minps(float,float);
int __builtin_ia32_maxps(float,float);
int __builtin_ia32_andps(float,float);
int __builtin_ia32_andnps(float,float);
int __builtin_ia32_orps(float,float);
int __builtin_ia32_xorps(float,float);
int __builtin_ia32_cmpeqss(float,float);
int __builtin_ia32_cmpltss(float,float);
int __builtin_ia32_cmpless(float,float);
int __builtin_ia32_cmpltss(float,float);
int __builtin_ia32_movss(float,float);
int __builtin_ia32_cmpless(float,float);
int __builtin_ia32_cmpneqss(float,float);
int __builtin_ia32_cmpnltss(float,float);
int __builtin_ia32_cmpnless(float,float);
int __builtin_ia32_cmpordss(float,float);
int __builtin_ia32_cmpunordss(float,float);
int __builtin_ia32_cmpeqps(float,float);
int __builtin_ia32_cmpltps(float,float);
int __builtin_ia32_cmpleps(float,float);
int __builtin_ia32_cmpgtps(float,float);
int __builtin_ia32_cmpgeps(float,float);
int __builtin_ia32_cmpneqps(float,float);
int __builtin_ia32_cmpnltps(float,float);
int __builtin_ia32_cmpnleps(float,float);
int __builtin_ia32_cmpngtps(float,float);
int __builtin_ia32_cmpngeps(float,float);
int __builtin_ia32_cmpordps(float,float);
int __builtin_ia32_cmpunordps(float,float);
int __builtin_ia32_comieq(float,float);
int __builtin_ia32_comilt(float,float);
int __builtin_ia32_comile(float,float);
int __builtin_ia32_comigt(float,float);
int __builtin_ia32_comige(float,float);
int __builtin_ia32_comineq(float,float);
int __builtin_ia32_ucomieq(float,float);
int __builtin_ia32_ucomilt(float,float);
int __builtin_ia32_ucomile(float,float);
int __builtin_ia32_ucomigt(float,float);
int __builtin_ia32_ucomige(float,float);
int __builtin_ia32_ucomineq(float,float);
int __builtin_ia32_cvtss2si(float);
int __builtin_ia32_cvtss2si64(float);
int __builtin_ia32_cvtps2pi(float);
int __builtin_ia32_cvttss2si(float);
int __builtin_ia32_cvttss2si64(float);
int __builtin_ia32_cvttps2pi(float);
int __builtin_ia32_cvtsi2ss(float,float);
int __builtin_ia32_cvtsi642ss(float,float);
int __builtin_ia32_cvtsi642ss(float,float);
int __builtin_ia32_cvtpi2ps(float,float);
int __builtin_ia32_cvtpi2ps(float,float);
int __builtin_ia32_movlhps(float,float);
int __builtin_ia32_cvtpi2ps(float,float);
int __builtin_ia32_movlhps(float,float);
int __builtin_ia32_cvtpi2ps(float,float);
int __builtin_ia32_movlhps(float,float);
int __builtin_ia32_movhlps(float,float);
int __builtin_ia32_cvtps2pi(float);
int __builtin_ia32_unpckhps(float,float);
int __builtin_ia32_unpcklps(float,float);
int __builtin_ia32_loadhps(float,int*);
int __builtin_ia32_storehps(int*,float);
int __builtin_ia32_movhlps(float,float);
int __builtin_ia32_movlhps(float,float);
int __builtin_ia32_loadlps(float,int*);
int __builtin_ia32_storelps(int*,float);
int __builtin_ia32_movmskps(float);
int __builtin_ia32_stmxcsr();
int __builtin_ia32_ldmxcsr(int);
int __builtin_ia32_loadups(float const*);
int __builtin_ia32_shufps(short,short,int);
int __builtin_ia32_vec_ext_v4sf(float,float);
int __builtin_ia32_vec_ext_v4sf(float,float);
int __builtin_ia32_storeups(float*,short);
int __builtin_ia32_pmaxsw(float,float);
int __builtin_ia32_pmaxub(float,float);
int __builtin_ia32_pminsw(float,float);
int __builtin_ia32_pminub(float,float);
int __builtin_ia32_pmovmskb(char);
int __builtin_ia32_pmulhuw(float,float);
int __builtin_ia32_maskmovq(char,char,char*);
int __builtin_ia32_pavgb(float,float);
int __builtin_ia32_pavgw(float,float);
int __builtin_ia32_psadbw(float,float);
int __builtin_ia32_movntq(unsigned long long*,unsigned long long);
int __builtin_ia32_movntps(float*,float);
int __builtin_ia32_sfence();


/* From: emmintrin.h */
int __builtin_ia32_movsd(double,double);
int __builtin_ia32_loadupd(double const *);
int __builtin_ia32_shufpd(double,double,int);
int __builtin_ia32_storeupd(double*,double);
int __builtin_ia32_vec_ext_v2df(double,int);
int __builtin_ia32_shufpd();
int __builtin_ia32_vec_ext_v4si(int,int);
int __builtin_ia32_vec_ext_v2di(long long,int);
int __builtin_ia32_addpd(double,double);
int __builtin_ia32_addsd(double,double);
int __builtin_ia32_subpd(double,double);
int __builtin_ia32_subsd(double,double);
int __builtin_ia32_mulpd(double,double);
int __builtin_ia32_mulsd(double,double);
int __builtin_ia32_divpd(double,double);
int __builtin_ia32_divsd(double,double);
int __builtin_ia32_sqrtpd(double);
int __builtin_ia32_sqrtsd(double);
int __builtin_ia32_minpd(double,double);
int __builtin_ia32_minsd(double,double);
int __builtin_ia32_maxpd(double,double);
int __builtin_ia32_maxsd(double,double);
int __builtin_ia32_andpd(double,double);
int __builtin_ia32_andnpd(double,double);
int __builtin_ia32_orpd(double,double);
int __builtin_ia32_xorpd(double,double);
int __builtin_ia32_cmpeqpd(double,double);
int __builtin_ia32_cmpltpd(double,double);
int __builtin_ia32_cmplepd(double,double);
int __builtin_ia32_cmpgtpd(double,double);
int __builtin_ia32_cmpgepd(double,double);
int __builtin_ia32_cmpneqpd(double,double);
int __builtin_ia32_cmpnltpd(double,double);
int __builtin_ia32_cmpnlepd(double,double);
int __builtin_ia32_cmpngtpd(double,double);
int __builtin_ia32_cmpngepd(double,double);
int __builtin_ia32_cmpordpd(double,double);
int __builtin_ia32_cmpunordpd(double,double);
int __builtin_ia32_cmpeqsd(double,double);
int __builtin_ia32_cmpltsd(double,double);
int __builtin_ia32_cmplesd(double,double);
int __builtin_ia32_cmpltsd(double,double);
int __builtin_ia32_cmplesd(double,double);
int __builtin_ia32_cmpneqsd(double,double);
int __builtin_ia32_cmpnltsd(double,double);
int __builtin_ia32_cmpnlesd(double,double);
int __builtin_ia32_cmpordsd(double,double);
int __builtin_ia32_cmpunordsd(double,double);
int __builtin_ia32_comisdeq(double,double);
int __builtin_ia32_comisdlt(double,double);
int __builtin_ia32_comisdle(double,double);
int __builtin_ia32_comisdgt(double,double);
int __builtin_ia32_comisdge(double,double);
int __builtin_ia32_comisdg(double,double);
int __builtin_ia32_comisdneq(double,double);
int __builtin_ia32_ucomisdeq(double,double);
int __builtin_ia32_ucomisdlt(double,double);
int __builtin_ia32_ucomisdle(double,double);
int __builtin_ia32_ucomisdgt(double,double);
int __builtin_ia32_ucomisdge(double,double);
int __builtin_ia32_ucomisdneq(double,double);
int __builtin_ia32_loaddqu(char const*);
int __builtin_ia32_storedqu(char*, char);
int __builtin_ia32_cvtdq2pd(int);
int __builtin_ia32_cvtdq2ps(int);
int __builtin_ia32_cvtpd2dq(double);
int __builtin_ia32_cvtpd2pi(double);
int __builtin_ia32_cvtpd2ps(double);
int __builtin_ia32_cvttpd2dq(double);
int __builtin_ia32_cvttpd2pi(double);
int __builtin_ia32_cvtpi2pd(int);
int __builtin_ia32_cvtps2dq(double);
int __builtin_ia32_cvttps2dq(double);
int __builtin_ia32_cvtps2pd(float);
int __builtin_ia32_cvtsd2si(double);
int __builtin_ia32_cvtsd2si64(double);
int __builtin_ia32_cvtsd2si64(double);
int __builtin_ia32_cvttsd2si(double);
int __builtin_ia32_cvttsd2si64(double);
int __builtin_ia32_cvtsd2ss(int,int);
int __builtin_ia32_cvtsi2sd(int,int);
int __builtin_ia32_cvtsi642sd(int,int);
int __builtin_ia32_cvtsi642sd(int,int);
int __builtin_ia32_cvtss2sd(int,int);
int __builtin_ia32_unpcklpd(int,int);
int __builtin_ia32_unpckhpd(double,double);
int __builtin_ia32_loadhpd(double, double const *);
int __builtin_ia32_loadlpd(double, double const *);
int __builtin_ia32_movmskpd(double);
int __builtin_ia32_packsswb128(short,short);
int __builtin_ia32_packssdw128(short,short);
int __builtin_ia32_packuswb128(short,short);
int __builtin_ia32_punpckhbw128(short,short);
int __builtin_ia32_punpckhwd128(int,int);
int __builtin_ia32_punpckhdq128(int,int);
int __builtin_ia32_punpckhqdq128(int,int);
int __builtin_ia32_punpcklbw128(int,int);
int __builtin_ia32_punpcklwd128(int,int);
int __builtin_ia32_punpckldq128(int,int);
int __builtin_ia32_punpcklqdq128(int,int);
int __builtin_ia32_paddb128(long long,long long);
int __builtin_ia32_paddw128(short,short);
int __builtin_ia32_paddd128(int,int);
int __builtin_ia32_paddq128(long long,long long);
int __builtin_ia32_paddsb128(long long,long long);
int __builtin_ia32_paddsw128(short,short);
int __builtin_ia32_paddusb128(char,char);
int __builtin_ia32_paddusw128(short,short);
int __builtin_ia32_psubb128(char,char);
int __builtin_ia32_psubw128(short,short);
int __builtin_ia32_psubd128(int,int);
int __builtin_ia32_psubq128(double,double);
int __builtin_ia32_psubsb128(char,char);
int __builtin_ia32_psubsw128(short,short);
int __builtin_ia32_psubusb128(char,char);
int __builtin_ia32_psubusw128(short,short);
int __builtin_ia32_pmaddwd128(short,short);
int __builtin_ia32_pmulhw128(short,short);
int __builtin_ia32_pmullw128(short,short);
int __builtin_ia32_pmuludq(int,int);
int __builtin_ia32_pmuludq128(int,int);
int __builtin_ia32_psllwi128(short,short);
int __builtin_ia32_pslldi128(int,int);
int __builtin_ia32_psllqi128(int,int);
int __builtin_ia32_psrawi128(short,short);
int __builtin_ia32_psradi128(short,short);
int __builtin_ia32_psrlwi128(short,short);
int __builtin_ia32_psrldi128(short,short);
int __builtin_ia32_psrlqi128(short,short);
int __builtin_ia32_psllw128(short,short);
int __builtin_ia32_pslld128(short,short);
int __builtin_ia32_psllq128(short,short);
int __builtin_ia32_psraw128(short,short);
int __builtin_ia32_psrad128(short,short);
int __builtin_ia32_psrlw128(short,short);
int __builtin_ia32_psrld128(short,short);
int __builtin_ia32_psrlq128(short,short);
int __builtin_ia32_pand128(int,int);
int __builtin_ia32_pandn128(int,int);
int __builtin_ia32_por128(int,int);
int __builtin_ia32_pxor128(int,int);
int __builtin_ia32_pcmpeqb128(char,char);
int __builtin_ia32_pcmpeqw128(short,short);
int __builtin_ia32_pcmpeqd128(short,short);
int __builtin_ia32_pcmpgtb128(char,char);
int __builtin_ia32_pcmpgtw128(short,short);
int __builtin_ia32_pcmpgtd128(int,int);
int __builtin_ia32_pcmpgtb128(char,char);
int __builtin_ia32_pcmpgtw128(short,short);
int __builtin_ia32_pmaxsw128(short,short);
int __builtin_ia32_pmaxub128(char,char);
int __builtin_ia32_pminsw128(short,short);
int __builtin_ia32_pminub128(char,char);
int __builtin_ia32_pmovmskb128(char);
int __builtin_ia32_pmulhuw128(short,short);
int __builtin_ia32_maskmovdqu(char,char,char*);
int __builtin_ia32_pavgb128(char,char);
int __builtin_ia32_pavgw128(short,short);
int __builtin_ia32_psadbw128(char,char);
int __builtin_ia32_movnti(int*,int);
int __builtin_ia32_movntdq(long long*,long long);
int __builtin_ia32_movntpd(double*,double);
int __builtin_ia32_clflush(void const *);
int __builtin_ia32_lfence();
int __builtin_ia32_mfence();

// Builtin functions specific to GNU 4.4.1 (likely 4.4.x)
int __builtin_ia32_psllwi(short,int);
int __builtin_ia32_pslldi(int,int);
int __builtin_ia32_psllqi(long long,int);
int __builtin_ia32_psrawi(short,int);
int __builtin_ia32_psradi(int,int);
int __builtin_ia32_psrlwi(int,int);
int __builtin_ia32_psrldi(int,int);
int __builtin_ia32_psrlqi(long long,int);

/* matching else for SKIP_ROSE_BUILTIN_DECLARATIONS */
#else

/* When compiling using -DSKIP_ROSE_BUILTIN_DECLARATIONS we need to have a variable
   defined in this faile so that we can locate the file in the AST and obtain the
   absolute path of this front-end specific header file.  This allows us to mark
   all IR nodes in the AST as being front-end specific if they originate from this
   header file. The following variable guarentees that a variable declaration from
   this file will exist when -DSKIP_ROSE_BUILTIN_DECLARATIONS is used.  The 
   -DSKIP_ROSE_BUILTIN_DECLARATIONS option is typically used to reduce the size of
   the AST and permit visualization of the whole AST using DOT (which can't layout
   a graph containing too many IR nodes).
 */
int __frontend_specific_variable_to_provide_header_file_path;

/* matching endif for SKIP_ROSE_BUILTIN_DECLARATIONS */
#endif

//I think these definitions will break on 32-bit systems
//The built-in declarations from gcc actually specify the bit length explictly
//e.g. int32_t __builtin_bswap32 (int32_t x)
//However, integer types such as int32_t are not defined here.
int __builtin_bswap32 (int x);
long int __builtin_bswap64 (long int x);

#endif /* !ROSE_USE_NEW_EDG_INTERFACE */
