#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <string.h>

#undef MALLOC
#define	MALLOC(name,len,type) {	\
	 (name) = (type *)malloc( (len) * sizeof(type) ); \
}

#undef CALLOC
#define	CALLOC(name,len,type) { \
	 (name) = (type *)calloc( (len) , sizeof(type) ); \
}

#undef REALLOC
#define	REALLOC(name,len,type) { \
	(name) = (type *)realloc((name),(len)*sizeof(type)); \
}

#undef FREE
#define	FREE(name) { \
	if ( (name) != NULL ) free( (name) ); \
	(name) = NULL; \
}
