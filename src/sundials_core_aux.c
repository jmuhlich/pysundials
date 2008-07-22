#include <stdio.h>
#include "sundials/sundials_types.h"
#include "sundials/sundials_smalldense.h"

int getsizeofrealtype() {
	return sizeof(realtype);
}

//int copyarray(realtype *dst, realtype *src, int size) {
//	memmove((void *)dst, (void *)src, size);
//}
