#include <cstdio>
#include "sun.h"
#include "coordutils.h"

main()
{
	char str[100], str2[100];
	double d;

	while(!feof(stdin))
	{
		scanf("%s", str);
		d = text2double(str);
		printf(" -> %f\n", d);
		double2text(d, str2);
		printf(" --> %s\n", str2);
	}
}
