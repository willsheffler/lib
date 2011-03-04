#include <stdio.h>
#include <stdlib.h>

int main (int argc, char const *argv[])
{
	float a = 3.1415;
	float b = 3.14159;

	int ia = *((int*)&a);
	int ib = *((int*)&b);
	
	int THRESH = 5; // depending on how close you want your floats to be...
	if( abs(ia-ib) <= THRESH ) {
		printf("close enough!\n");
	}
	
	printf("result: %f %f %d\n",a,b,ia-ib);
	
	return 0;
}
