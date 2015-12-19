#ifndef CONFIG_H_
#define CONFIG_H_

#ifdef _WIN_X64
	#pragma warning(disable:4267)
#endif
#pragma warning(disable:4996)

#define USE_DOUBLE_AS_FLOAT

#ifndef USE_DOUBLE_AS_FLOAT
	typedef float  piFloat;
	#define  piFloatEpsilon		1e-3
#else
	typedef double piFloat;
	#define  piFloatEpsilon		1e-5
#endif

//#define piFloatEpsilon std::numeric_limits<piFloat>::epsilon()

#define DEFAULT_BUFFER_PIXEL_THRESHOLD	10

#endif
