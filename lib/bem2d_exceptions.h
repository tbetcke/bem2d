
#ifndef _EXCEPTIONS_H
#define	_EXCEPTIONS_H

// bem2d exception classes

namespace bem2d
{

class ArrayMismatch
{
};

class LapackError
{
};

class IntervalException
{
};

class ParameterException
{

};

class SizeError {};

#ifdef BEM2DMPI
	class ScaLapackError
	{
	};
#endif
	
}



#endif	/* _EXCEPTIONS_H */

