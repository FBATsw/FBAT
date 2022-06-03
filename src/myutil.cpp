#include <stdio.h>
#include "myutil.h"

char* GetWord(char *s, int whichword, char *w) {
	int i;
	char *c;
	
	while (*s==' ')
		s++;
	for (i=1; *s && i<whichword; s++)
		if ((*s==' ' && *(s-1)!=' ') || (*s=='\t' && *(s-1)!=' '))
			i++;

	if (i<whichword)
		return 0;
	
	if (*(s-1)==' ') {
		while (*s==' ')
			s++;
		if (*s=='\t')
			s++;
	}
	
	if (w==0)
		return s;
		
	c = s;
	
	for (i=0; *s && *s!=' ' && *s!='\t' && *s!='\n' && *s!='\r'; i++)
		*w++ = *s++;
	*w = 0;
	
	return c;

}

int WordCount(char *s) {
	
	int cnt, t;
	
	for (t=0; *s==' ' || *s=='\t'; s++)
		if (*s=='\t')
			t++;
	cnt = t;
	while (*s && *s!='\n' && *s!='\r') {
		while (*s && *s!='\n' && *s!='\r' && *s!=' ' && *s!='\t')
			s++;
		for (t=0; *s==' ' || *s=='\t'; s++)
			if (*s=='\t')
				t++;
		cnt += t?t:1;
	}
	
	return t?cnt+1:cnt;

}

// print out p values, if p<10^-5, use scientific notation
void pval2str(double p, char *s) {
	
	if (p<0.00001)
		sprintf(s, "%.2le",p);
	else
		sprintf(s, "%.6lf",p);

}



	