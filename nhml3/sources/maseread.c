#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/*
To read a mase file

call it as follows:
	const int maxnseqs=100;
	char *seq[maxnseqs], *seqname[maxnseqs], masefname[40];
	char *comments[maxnseqs] or char **comments=NULL;
	int totseqs;
	totseqs=readmaseseqs(masefname,seq,seqname,comments,maxnseqs);
	if(totseqs==0) <nothing read>;

masefname (input): name of mase file to read
seq (input/output): returned loaded with array of sequences read 
	(regions only are loaded when they are given in mase file
seqname (input/output): returned with array of sequence names
comments (input/output): if transmitted NULL, unused
          else, loaded with comments of each sequence (with their \n)
maxnseqs (input): max number of sequences = size of arrays *seq and *seqname
returned value: # of sequences found in file, or 0 if error

If regions are not given in mase file, the max length of sequences
is MAXLENSEQ below. If regions appear, memory is dynamically allocated
*/

int readmaseseqs(char *masefname, char **seq, char **seqname, char ** comments,
							 int maxnseqs)
{
/* long max des seqs sans regions; avec regions la longueur s'ajuste seule */
#define MAXLENSEQ 10000
#define MAXLENCOM 5000 /* long max des commentaires */
int lline=140;
FILE *masef;
char line[140], *i, rep[10], *base;
int nreg, treg, *endpoints, l, lenseqs=0, minreg, lpre, lseq, l2, totseqs=-1,
	 curreg;
if( (masef=fopen(masefname,"r")) == NULL) {
	fprintf(stderr,"File not found:%s\n",masefname);
	return 0;
	}
endpoints=NULL;
treg=0;
do	{
	if(fgets(line,lline,masef)==NULL)goto fini;
	if(treg==0 && strncmp(line,";;#",3)==0) {
		printf("\n%sUse this set of regions? ([y]/n) ",line);
		*rep='\0'; gets(rep);
		if(*rep=='n' || *rep=='N') *(line+2)='\0';
		}
	if(treg==0 && strncmp(line,";;#",3)==0) {
		i=strchr(line,'=');
		sscanf(i+1,"%d",&treg);
		endpoints= (int *) malloc(2*treg*sizeof(int));
		if(endpoints==NULL)goto nomem;
		nreg=0;
		while(nreg<=2*treg-1) {
			char separ[5], *j;
			fgets(line,lline,masef);
			strcpy(separ,";, \n"); j=line;
			while( (i=strtok(j,separ)) != NULL) {
				sscanf(i,"%d",&l);
				if(nreg>0 && l<=*(endpoints+nreg-1)) {
					fprintf(stderr,
			"Region endpoints are not in increasing order: %d\n",l);
					goto fini;
					}
				*(endpoints+nreg++)=l-1;
				j=NULL;
				}
			}
		}
	}
while(strncmp(line,";;",2)==0);
if(endpoints!=NULL) {
	printf("\nEndpoints of regions used:\n");
	for (l=1;l<2*treg;l+=2) {
		printf(" %d,%d",*(endpoints+l-1)+1,*(endpoints+l)+1);
		if( (l+1)%12 == 0)printf("\n");
		lenseqs += *(endpoints+l)-*(endpoints+l-1)+1;
		}
	printf("\n");
	}
else	{
	printf("\nComplete sequences are used.\n");
	lenseqs=MAXLENSEQ;
	treg=1; 
	if( (endpoints=(int *)malloc(2*sizeof(int)) )==NULL)goto nomem;
	*endpoints=0; *(endpoints+1)=lenseqs-1;
	}

i=line;
while(i!=NULL){
	if(totseqs>=maxnseqs-1) {
		printf("Reading of mase file stopped after\
 %d sequences\n",totseqs+1); 
		goto fini;
		}
	totseqs++;
	if(comments!=NULL) {
		if( (comments[totseqs]=(char *)malloc(MAXLENCOM+1)) ==
							 NULL)goto nomem;
		strcpy(comments[totseqs],line);
		lpre=strlen(line); l=MAXLENCOM;
		while(*fgets(line,lline,masef)==';') {
			lseq=strlen(line);
			if(lpre+lseq <= l) {
				strcpy(comments[totseqs]+lpre,line);
				lpre += lseq;
				}
			else l=lpre-1;
			}
		if(lpre<MAXLENCOM)
		   comments[totseqs]=(char *)realloc(comments[totseqs],lpre+1);
		}
	else	while(*fgets(line,lline,masef)==';');
	if( (seq[totseqs]=(char *)malloc(lenseqs+1)) ==NULL)goto nomem;
	l=strlen(strtok(line," \n"));
	if( (seqname[totseqs]=(char *)malloc(l+1)) == NULL)goto nomem;
	strcpy(seqname[totseqs],line);
	minreg=*endpoints; /* start of region that remains to be put in seq */
	curreg=0; /* current region processed: 0, 2, 4... */
	lpre=-1; /* what was read in previous lines */
	lseq=-1; /* what is already put in seq */
	while( (i=fgets(line,lline,masef))!= NULL && *i != ';' ) {
	   line[strlen(line)-1]='\0';
	   l2=strlen(line);
	   while(l2>0 && line[l2-1]==' ')l2--;
	   while(curreg<=2*treg-1 && minreg<=l2+lpre) {
		l=*(endpoints+curreg+1);
		if(l>l2+lpre) l=l2+lpre;
		strncpy(seq[totseqs]+lseq+1,line+minreg-lpre-1,l-minreg+1);
		lseq+=l-minreg+1;
		if(l==*(endpoints+curreg+1)) {
			curreg+=2;
			if(curreg<=2*treg-1)minreg=*(endpoints+curreg);
			}
		else	{
			minreg=l+1;
			}
		}
	   lpre+=l2;
	   }
	seq[totseqs][++lseq]='\0';
	if(lseq<lenseqs)seq[totseqs]=(char *)realloc(seq[totseqs],lseq+1);
/* convert all to upper case and change U to T */
	base=seq[totseqs];
	while (*base != '\0') {
		if(isalpha(*base)) {
			if(islower(*base)) *base=toupper(*base);
			if(*base=='U') *base='T';
			}
		else
			*base='-';
		base++;
		}
	}
fini:
fclose(masef);
return totseqs+1;
nomem:
fprintf(stderr,"Not enough memory!\n");
goto fini;
}


int
readmasemuet(char *masefname, char **seq, char **seqname, char **comments,
	     int maxnseqs)
{
	/*
	 * long max des seqs sans regions; avec regions la longueur s'ajuste
	 * seule
	 */
#define MAXLENSEQ 10000
#define MAXLENCOM 5000		/* long max des commentaires */
#define lline 140       
	FILE           *masef;
	char            line[lline], *i, rep[10], *base;
	int             nreg, treg, *endpoints, l, lenseqs = 0, minreg, lpre,
	                lseq, l2, totseqs = -1, curreg;
	if ((masef = fopen(masefname, "r")) == NULL) {
		fprintf(stderr, "File not found:%s\n", masefname);
		return 0;
	}
	endpoints = NULL;
	treg = 0;
	do {
		if (fgets(line, lline, masef) == NULL)
			goto fini;
	}
	while (strncmp(line, ";;", 2) == 0);
	if (endpoints != NULL) {
		for (l = 1; l < 2 * treg; l += 2) {
			printf(" %d,%d", *(endpoints + l - 1) + 1, *(endpoints + l) + 1);
			if ((l + 1) % 12 == 0)
				printf("\n");
			lenseqs += *(endpoints + l) - *(endpoints + l - 1) + 1;
		}
	} else {
		lenseqs = MAXLENSEQ;
		treg = 1;
		if ((endpoints = (int *) malloc(2 * sizeof(int))) == NULL)
			goto nomem;
		*endpoints = 0;
		*(endpoints + 1) = lenseqs - 1;
	}

	i = line;
	while (i != NULL) {
		if (totseqs >= maxnseqs - 1) {
			printf("Reading of mase file stopped after\
 %d sequences\n", totseqs + 1);
			goto fini;
		}
		totseqs++;
		if (comments != NULL) {
			if ((comments[totseqs] = (char *) malloc(MAXLENCOM + 1)) ==
			    NULL)
				goto nomem;
			strcpy(comments[totseqs], line);
			lpre = strlen(line);
			l = MAXLENCOM;
			while (*fgets(line, lline, masef) == ';') {
				lseq = strlen(line);
				if (lpre + lseq <= l) {
					strcpy(comments[totseqs] + lpre, line);
					lpre += lseq;
				} else
					l = lpre - 1;
			}
			if (lpre < MAXLENCOM)
				comments[totseqs] = (char *) realloc(comments[totseqs], lpre + 1);
		} else
			while (*fgets(line, lline, masef) == ';');
		if ((seq[totseqs] = (char *) malloc(lenseqs + 1)) == NULL)
			goto nomem;
		l = strlen(strtok(line, " \n"));
		if ((seqname[totseqs] = (char *) malloc(l + 1)) == NULL)
			goto nomem;
		strcpy(seqname[totseqs], line);
		minreg = *endpoints;	/* start of region that remains to be
					 * put in seq */
		curreg = 0;	/* current region processed: 0, 2, 4... */
		lpre = -1;	/* what was read in previous lines */
		lseq = -1;	/* what is already put in seq */
		while ((i = fgets(line, lline, masef)) != NULL && *i != ';') {
			line[strlen(line) - 1] = '\0';
			l2 = strlen(line);
			while (l2 > 0 && line[l2 - 1] == ' ')
				l2--;
			while (curreg <= 2 * treg - 1 && minreg <= l2 + lpre) {
				l = *(endpoints + curreg + 1);
				if (l > l2 + lpre)
					l = l2 + lpre;
				strncpy(seq[totseqs] + lseq + 1, line + minreg - lpre - 1, l - minreg + 1);
				lseq += l - minreg + 1;
				if (l == *(endpoints + curreg + 1)) {
					curreg += 2;
					if (curreg <= 2 * treg - 1)
						minreg = *(endpoints + curreg);
				} else {
					minreg = l + 1;
				}
			}
			lpre += l2;
		}
		seq[totseqs][++lseq] = '\0';
		if (lseq < lenseqs)
			seq[totseqs] = (char *) realloc(seq[totseqs], lseq + 1);
		/* convert all to upper case and change U to T */
		base = seq[totseqs];
		while (*base != '\0') {
			if (isalpha(*base)) {
				if (islower(*base))
					*base = toupper(*base);
				if (*base == 'U')
					*base = 'T';
			} else
				*base = '-';
			base++;
		}
	}
fini:
	fclose(masef);
	return totseqs + 1;
nomem:
	fprintf(stderr, "Not enough memory!\n");
	goto fini;
}
