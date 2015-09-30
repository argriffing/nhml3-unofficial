#include "nhmlg.h"
#define MAXNTREE 100

tree* s_tree;
compute_option compute_o;
noeud root;
int gamma_nbcl;

void printmemory_init();
void printmemory(char []);
void init_cas();
void *check_alloc(int nbrelt, int sizelt);
void init_bl(tree* stree);
void init4bases_constant(tree* s_tree);
void init4bases_simplemean(noeud from, noeud nd, int lgseq, double weightsum, double* weight);
void set_listgc(tree* s_tree, noeud nd);
void reset_like(tree* s_tree);
void boundaries(tree* s_tree, int* bounded);
void reset_tree(tree* s_tree);
double minus_log_likelihood(double *paramprov);
void derivatives(tree* s_tree, int cross);
double simplex(double* initpoint, double bound, double* finalvect);
int invmat(double** mat, int n, double** invmat);
int mindepth (noeud comesfrom, noeud node);
void makelistbr_unrooted(noeud from, noeud nd, branche* br, int* alivebr,int nbbranch);
void set_var_param(tree* s_tree, compute_option comp);
double alsurbet_kim(char* seq1, char* seq2, double* freq);
void organize_tree(noeud from, noeud nd);
void organize_listbr(tree* s_tree);
void setpattern(tree* s_tree);
int setunderlyingbr(noeud nd, branche br, int i);
int setunderlyinggc(noeud nd, int nb);
void set_numerobr(noeud nd, branche* listbr, int nbbr);
void set_node_deduced_comp(noeud root);
noeud* ctos(char* ctree, int* notu, int *tree_is_rooted);
void stoc(noeud* arbre_s, int racine, int notu, char* ctree, int nodegc);
void refresh(char** seq, int nbseq, int option);
double gg95(char* seq1, char* seq2, double alpha);
double compute(tree* s_tree, init_option init_o, compute_option compute_opt, converge_option converge_o, print_option print_o);
char* readtree(FILE* treef, int nbsp);
void getoptions(options* opt_arg, FILE* in_opt);




void restore_lengths(tree* s_tree){
 
  int i, nbseq;
  noeud nd1, nd2;
  branche br;
  double l;
  
  nbseq=s_tree->nbseq;

  for(i=0;i<2*nbseq-3;i++){
    br=s_tree->listbr[i];
    l=s_tree->lgbr_init[i];
    br->length=l;
    nd1=br->bout1;
    nd2=br->bout2;

    if(nd2==nd1->v3){  
      nd1->l3=l; 
      if(nd1==nd2->v1) nd2->l1=l;
      else nd2->l2=l;
    }
    else if(nd1==nd2->v3){
      nd2->l3=l; 
      if(nd2==nd1->v1) nd1->l1=l;
      else nd1->l2=l;
    }
    else if(nd1->v3==root && nd2->v3==root){
      nd1->l3=l; nd2->l3=l; 
      root->l1=l*s_tree->froot; root->l2=l*(1.-s_tree->froot);
    }
    else {
      printf("unexpected unrooted tree\n");
      exit(EXIT_FAILURE);
    }
  }


}





void free_node(noeud nd, int lgseq, int nbparam, int leaf, int nbcl){

  int j;

  free(nd->underlyingbr);
  free(nd->underlyinggc);
  free(nd->seq);
  if(!leaf) free(nd->patternson1);
  if(!leaf) free(nd->patternson2);

  for(j=0;j<nbcl;j++) free(nd->x[j]);
  free(nd->x);
  for(j=0;j<nbcl;j++) free(nd->dx[j]);
  free(nd->dx);
  for(j=0;j<nbcl;j++) free(nd->d2x[j]);
  free(nd->d2x);
  for(j=0;j<nbcl;j++) free(nd->p1[j]);
  free(nd->p1);
  for(j=0;j<nbcl;j++) free(nd->p2[j]);
  free(nd->p2);
  for(j=0;j<nbcl;j++) free(nd->dp1[j]);
  free(nd->dp1);
  for(j=0;j<nbcl;j++) free(nd->dp2[j]);
  free(nd->dp2);
  for(j=0;j<nbcl;j++) free(nd->d2p1[j]);
  free(nd->d2p1);
  for(j=0;j<nbcl;j++) free(nd->d2p2[j]);
  free(nd->d2p2);
  for(j=0;j<nbcl;j++) free(nd->rr1[j]);
  free(nd->rr1);
  for(j=0;j<nbcl;j++) free(nd->rr2[j]);
  free(nd->rr2);
  for(j=0;j<nbcl;j++) free(nd->drr1[j]);
  free(nd->drr1);
  for(j=0;j<nbcl;j++) free(nd->drr2[j]);
  free(nd->drr2);
  for(j=0;j<nbcl;j++) free(nd->d2rr1[j]);
  free(nd->d2rr1);
  for(j=0;j<nbcl;j++) free(nd->d2rr2[j]);
  free(nd->d2rr2);
  free(nd);

}



void free_tree(tree* s_tree, int nbcl){

  int i, j, nbseq, nbparam, lgseq;

  nbseq=s_tree->nbseq;
  lgseq=s_tree->lgseq;
  nbparam=s_tree->nb_var_param;

  for(i=0;i<2*nbseq-2;i++)
    free_node(s_tree->node[i], lgseq, nbparam, (i<nbseq), nbcl);
  for(i=0;i<2*nbseq-3;i++)
    free(s_tree->listbr[i]);
  free(s_tree->listbr);
  free(s_tree->listgc);
  free(s_tree->lgbr_init);
  free(s_tree->deriv2);
  free(s_tree->deriv);
  for(i=0;i<nbcl;i++) free(s_tree->classlike_s[i]);
  free(s_tree->classlike_s);
  free(s_tree->like_s);
  free(s_tree->likeasrv_s);
  for(j=0;j<nbcl;j++){
    for(i=0;i<nbparam;i++){
      free(s_tree->classdlike_s[j][i]);
      free(s_tree->classd2like_s[j][i]);
    }
    free(s_tree->classdlike_s[j]);
    free(s_tree->classd2like_s[j]);
  }
  free(s_tree->classdlike_s);
  free(s_tree->classd2like_s);
  for(i=0;i<nbparam;i++){
    free(s_tree->dlike_s[i]);
    free(s_tree->dlikeasrv_s[i]);
    free(s_tree->d2like_s[i]);
    free(s_tree->d2likeasrv_s[i]);
  }
  free(s_tree->alivebr);
  free(s_tree->dlike_s);
  free(s_tree->dlikeasrv_s);
  free(s_tree->d2like_s);
  free(s_tree->d2likeasrv_s);
  free(s_tree->weight);
  free(s_tree->var_param);
  free(s_tree->param_nature);
  for(i=0;i<nbseq;i++) free(s_tree->dist[i]);
  free(s_tree->dist);
  free(s_tree->sitetopatt);
  
}



void alloc_node(noeud nd, int lgseq, int nbparam, int nbseq, int leaf, int cross, int nbcl){

  int j, i;

  nd->underlyingbr=(int*)check_alloc(2*nbseq-3, sizeof(int));
  nd->underlyinggc=(int*)check_alloc(2*nbseq-2, sizeof(int));
  nd->seq=(char*)check_alloc(lgseq, sizeof(char));
  if(!leaf) nd->patternson1=(int*)check_alloc(lgseq, sizeof(int));
  if(!leaf) nd->patternson2=(int*)check_alloc(lgseq, sizeof(int));

  nd->x=(like**)check_alloc(nbcl, sizeof(like*));
  nd->dx=(like**)check_alloc(nbcl, sizeof(like*));    
  nd->d2x=(like**)check_alloc(nbcl, sizeof(like*));    
  for(j=0;j<nbcl;j++) {
    nd->x[j]=(like*)check_alloc(lgseq, sizeof(like));
    nd->dx[j]=(like*)check_alloc(lgseq, sizeof(like));
    nd->d2x[j]=(like*)check_alloc(lgseq, sizeof(like));
  }

/*
  nd->x2sites=check_alloc(nbcl, sizeof(double**));
  for(j=0;j<nbcl;j++){
    nd->x2sites[j]=check_alloc(4, sizeof(double*));
    for(i=0;i<4;i++) nd->x2sites[j][i]=check_alloc(4, sizeof(double));
  }
*/

  nd->p1=(proba**)check_alloc(nbcl, sizeof(proba*));
  nd->p2=(proba**)check_alloc(nbcl, sizeof(proba*));
  nd->dp1=(dproba**)check_alloc(nbcl, sizeof(dproba*));
  nd->dp2=(dproba**)check_alloc(nbcl, sizeof(dproba*));
  nd->d2p1=(dproba**)check_alloc(nbcl, sizeof(dproba*));
  nd->d2p2=(dproba**)check_alloc(nbcl, sizeof(dproba*));
  for(j=0;j<nbcl;j++){
    nd->p1[j]=(proba*)check_alloc(nbcl, sizeof(proba));
    nd->p2[j]=(proba*)check_alloc(nbcl, sizeof(proba));
    nd->dp1[j]=(dproba*)check_alloc(nbcl, sizeof(dproba));
    nd->dp2[j]=(dproba*)check_alloc(nbcl, sizeof(dproba));
    nd->d2p1[j]=(dproba*)check_alloc(nbcl, sizeof(dproba));
    nd->d2p2[j]=(dproba*)check_alloc(nbcl, sizeof(dproba));
  }

  nd->rr1=check_alloc(nbcl, sizeof(double*));
  nd->rr2=check_alloc(nbcl, sizeof(double*));
  nd->drr1=check_alloc(nbcl, sizeof(double*));
  nd->drr2=check_alloc(nbcl, sizeof(double*));
  nd->d2rr1=check_alloc(nbcl, sizeof(double*));
  nd->d2rr2=check_alloc(nbcl, sizeof(double*));

  for(i=0;i<nbcl;i++){
    nd->rr1[i]=check_alloc(nbcl, sizeof(double));
    nd->rr2[i]=check_alloc(nbcl, sizeof(double));
    nd->drr1[i]=check_alloc(nbcl, sizeof(double));
    nd->drr2[i]=check_alloc(nbcl, sizeof(double));
    nd->d2rr1[i]=check_alloc(nbcl, sizeof(double));
    nd->d2rr2[i]=check_alloc(nbcl, sizeof(double));
  }


}




double maxlike(int nbseq, char** seq, char** seqname, char* c_tree, options opt, char* ctree1, char* ctree2){

  int i, j, k, l, ii, nbseq2, lgseq, lgsite, *drapeau, *drapeau2, cont=0, tree_is_rooted, nbparam, nbvrai, nbcompute, bl;
  char **orderedseq, racine, *site;
  double *lgbi, *lgbp, *weight, lkh, maxlkh=-1.e100, prov, rem_titv, rem_root, *freq;
  static int nbtree=0;
  FILE* outfile;

/* printmemory("deb maxlike"); */

  gamma_nbcl=opt->init->GAMMA_NCL;

  /*  READ TREE, GET NB OF TAXA */

  s_tree=(tree*)check_alloc(1, sizeof(tree));
  lgseq=(int)strlen(seq[0]);
  s_tree->alivebr=(int*)check_alloc(2*nbseq, sizeof(int));
  s_tree->node=ctos(c_tree, &nbseq2, &tree_is_rooted);
  if(strchr(c_tree, ':')) bl=1; else bl=0;


  if (nbseq2>nbseq) {
    printf("More taxa in tree file than in sequence file\n");
    exit(EXIT_FAILURE);
  }


	/* ALLOCATE AND INIT S_TREE */

  s_tree->nbseq=nbseq=nbseq2;
  s_tree->seq=seq;
  s_tree->names=seqname;

  s_tree->listgc=(double*)check_alloc(2*nbseq, sizeof(double));
  s_tree->listbr=(branche*)check_alloc(2*nbseq, sizeof(branche));
  s_tree->lgbr_init=(double*)check_alloc(2*nbseq-3, sizeof(double));
  for(i=0;i<2*nbseq-2;i++){
    s_tree->listbr[i]=(struct branche*)check_alloc(1, sizeof(struct branche));
  }

  makelistbr_unrooted(NULL, s_tree->node[0], s_tree->listbr, s_tree->alivebr, 2*nbseq-3);
  for(i=0;i<2*nbseq-3;i++) s_tree->lgbr_init[i]=s_tree->listbr[i]->length;
  for(ii=0;ii<2*nbseq-3;ii++) s_tree->alivebr[ii]=1;
  s_tree->nbalivebr=0;
  for(i=0;i<2*nbseq-3;i++) s_tree->nbalivebr+=s_tree->alivebr[i];


  set_var_param(s_tree, opt->compute);

  nbparam=s_tree->nb_var_param;

  s_tree->classlike_s=(double**)check_alloc(gamma_nbcl, sizeof(double));
  for(i=0;i<gamma_nbcl;i++)
    s_tree->classlike_s[i]=(double*)check_alloc(lgseq, sizeof(double));
  s_tree->like_s=(double*)check_alloc(lgseq, sizeof(double));
  s_tree->likeasrv_s=(double*)check_alloc(lgseq, sizeof(double));
  s_tree->dlike_s=(double**)check_alloc(nbparam, sizeof(double*));
  s_tree->dlikeasrv_s=(double**)check_alloc(nbparam, sizeof(double*));
  for(i=0;i<nbparam;i++){
    s_tree->dlike_s[i]=(double*)check_alloc(lgseq, sizeof(double));
    s_tree->dlikeasrv_s[i]=(double*)check_alloc(lgseq, sizeof(double));
  }
  s_tree->d2like_s=(double**)check_alloc(nbparam, sizeof(double*));
  s_tree->d2likeasrv_s=(double**)check_alloc(nbparam, sizeof(double*));
  for(i=0;i<nbparam;i++){
    s_tree->d2like_s[i]=(double*)check_alloc(lgseq, sizeof(double));
    s_tree->d2likeasrv_s[i]=(double*)check_alloc(lgseq, sizeof(double));
  }

  s_tree->classdlike_s=(double***)check_alloc(gamma_nbcl, sizeof(double**));
  s_tree->classd2like_s=(double***)check_alloc(gamma_nbcl, sizeof(double**));
  for(i=0;i<gamma_nbcl;i++){ 
    s_tree->classdlike_s[i]=(double**)check_alloc(nbparam, sizeof(double*));
    s_tree->classd2like_s[i]=(double**)check_alloc(nbparam, sizeof(double*));
    for(j=0;j<nbparam;j++){
      s_tree->classdlike_s[i][j]=(double*)check_alloc(lgseq, sizeof(double));
      s_tree->classd2like_s[i][j]=(double*)check_alloc(lgseq, sizeof(double));
    }
  }
  
  s_tree->deriv=(double*)check_alloc(nbparam, sizeof(double));
  s_tree->deriv2=(double*)check_alloc(nbparam, sizeof(double));

/* printmemory("fin tree alloc"); */

  for(i=0;i<2*nbseq-2;i++)
    alloc_node(s_tree->node[i], lgseq, nbparam, nbseq, (i<nbseq), 0, gamma_nbcl);
  alloc_node(root, lgseq, nbparam, nbseq, 0, 0, gamma_nbcl);

/* printmemory("fin node alloc"); */



	/* INIT LIKELIHOOD AT TIPS */

  for(i=0;i<nbseq;i++){
    for(j=0;j<4;j++){ /* 4 patterns for tip nodes */
      for(k=0;k<4;k++){ /* 4 possible states */
        for(l=0;l<gamma_nbcl;l++){
	  if(j==k) s_tree->node[i]->x[l][j][k]=1.;
	  else s_tree->node[i]->x[l][j][k]=0.;
        }
      }
    }
  }

       /* CHECK SEQUENCE FILE */

  for(i=nbseq2;i<2*nbseq2-2;i++) sprintf(s_tree->node[i]->nom, "int%d", i-nbseq2+1);

  drapeau=(int*)check_alloc(nbseq, sizeof(int));
  drapeau2=(int*)check_alloc(nbseq, sizeof(double));
  for(i=0;i<nbseq2;i++) drapeau2[i]=0;
  orderedseq=check_alloc(MAXNSP, sizeof(char*));
  

  for(i=0;i<nbseq2;i++){
    for(j=0;j<nbseq;j++){
      if(!drapeau[j] && strncmp(s_tree->node[i]->nom, seqname[j], MAXLNAME)==0){
        drapeau[j]=drapeau2[i]=1;
	orderedseq[i]=seq[j];
        break;
      }
    }
    if(drapeau2[i]==0){
      printf("Taxon %s not found in sequence file.\n", s_tree->node[i]->nom);
      exit(EXIT_FAILURE);
    }
  }
  s_tree->nbseq=nbseq=nbseq2;


	/* SET "PATTERNS" (relationships between sites) */

  s_tree->sitetopatt=check_alloc(lgseq, sizeof(int));
  weight=(double*)check_alloc(lgseq, sizeof(double));
  for(i=0;i<lgseq;i++) weight[i]=0.;
  site=(char*)check_alloc(nbseq, sizeof(char)); 
  lgsite=0;
  for(i=0;i<lgseq;i++){
    cont=0;
    for(j=0;j<nbseq;j++) site[j]=orderedseq[j][i];
    for(k=0;k<lgsite;k++){
      if(samesite(site, s_tree, k)){
        weight[k]++;
	cont=1;
        break;
      }
    }
    if(cont){
      s_tree->sitetopatt[i]=k;
      continue; 
    }
    for(j=0;j<nbseq;j++)
      s_tree->node[j]->seq[lgsite]=site[j];
    weight[lgsite]=1.;
    s_tree->sitetopatt[i]=lgsite;
    lgsite++;
  }

  free(drapeau); free(drapeau2); free(site);

  s_tree->lgseq=lgsite;
  s_tree->weight=weight;
  s_tree->weightsum=lgseq;

  /* printf("%d species, %d sites, %d patterns\n", nbseq, lgseq, lgsite);*/


	/* AMMEND INPUT */
       	
		/* ts/tv */

  rem_titv=opt->init->INIT_TITV;
  if(opt->init->INIT_TITV<=0.){
    freq=(double*)check_alloc(16, sizeof(double));
    opt->init->INIT_TITV=0.;
    nbvrai=nbseq*(nbseq-1)/2;
    for(i=0;i<nbseq;i++){
      for(j=i+1;j<nbseq;j++){
	prov=alsurbet_kim(seq[i], seq[j], freq);
  	if(prov<-0.5) nbvrai--; else opt->init->INIT_TITV+=prov;
      }
    }
    free(freq);
    opt->init->INIT_TITV/=nbvrai;
  }

		/* root */

  rem_root=opt->init->INIT_ROOT;
  if(opt->init->INIT_ROOT<=0.){
    opt->init->INIT_ROOT=root->l1/(root->l1+root->l2);
    if(!bl) opt->init->INIT_ROOT=0.5;
  }

		/* branch lengths */

  s_tree->dist=(double**)check_alloc(nbseq, sizeof(double*));
  for(i=0;i<nbseq;i++)
    s_tree->dist[i]=(double*)check_alloc(nbseq, sizeof(double));

  if(strcmp(opt->init->INIT_LENGTH, "REDO")==0){
    for(i=0;i<nbseq;i++){
      for(j=i;j<nbseq;j++){
        if(i==j) s_tree->dist[i][j]=0.;
        else s_tree->dist[i][j]=s_tree->dist[j][i]=gg95(seq[i], seq[j], opt->init->INIT_TITV);
      }
    }
  }

	/* ROOT TREE AND ORGANIZE */

  organize_tree(NULL, root);
  organize_listbr(s_tree);
  setpattern(s_tree);
  for(ii=0;ii<2*nbseq-3;ii++)
    setunderlyingbr(root, s_tree->listbr[ii], ii);
  setunderlyinggc(root, s_tree->nbseq);
  set_numerobr(root, s_tree->listbr, 2*nbseq-3);

  for(i=0;i<nbseq;i++){
    s_tree->node[i]->alive1=0;
    s_tree->node[i]->alive2=0;
    s_tree->node[i]->alive3=1;
  }
  for(i=nbseq;i<2*nbseq-2;i++){
    s_tree->node[i]->alive1=1;
    s_tree->node[i]->alive2=1;
    s_tree->node[i]->alive3=1;
  }
  root->alive1=root->alive2=1; root->alive3=-1;
  for(ii=0;ii<2*nbseq-3;ii++) s_tree->alivebr[ii]=1;


	/* CALL COMPUTE */

  if(opt->init->NBRANDOM>0)
    nbcompute=opt->init->NBRANDOM;
  else 
    nbcompute=1;


  for(l=0;l<nbcompute;l++){
    lkh=compute(s_tree, opt->init, opt->compute, opt->converge, opt->print);
    if(lkh>maxlkh) maxlkh=lkh;
    if(!opt->print->PRINT1 && !opt->print->PRINT2) printf("%f\n", lkh);
  }


  if(ctree1)
    stoc(s_tree->node, nbseq, 1, ctree1, 0);

  if(ctree2)
    stoc(s_tree->node, nbseq, 1, ctree2, 1);

  if(opt->print->EVAL_OUT){
    outfile=fopen("detailed_out", "a");
    nbtree++;
    fprintf(outfile, "%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\n", nbtree, ctree1, maxlkh, s_tree->titv, s_tree->GCanc, s_tree->gamma_shp, s_tree->covar, s_tree->pi);
    fclose(outfile);
  }


  opt->init->INIT_TITV=rem_titv;
  opt->init->INIT_ROOT=rem_root;


  /* FREE */

  free_tree(s_tree, gamma_nbcl);
  free(s_tree);
  free_node(root, lgseq, nbparam, 0, gamma_nbcl);

  free(orderedseq);

  return maxlkh;

}



main(int argc, char** argv){

/* 		A->0 ; C->1 ; G->2 ; T->3 		*/

  int i, j, nbseq, nbseqvrai, nbtree, nbdataset, maxtreesize, allcouples=0;
  char nomfinseq[100], nomfintree[100], nomfopt[100], *seq[MAXNSP*MAXDATASET], *seqname[MAXNSP*MAXDATASET], *comments[MAXNSP*MAXDATASET], muet, print1;
  char* alltrees, **c_tree, *prov, *prov2, *ctree1, *ctree2;
  double **maxl;
  FILE *treefile, *optfile, *outfile1, *outfile2, *in;
  FILE* dutheil;
  options opt;

/* printmemory_init(); */

  init_cas();

  srand48(seed);


   /** input sequences **/


  if(argc<2){
    muet=0;
    while(1){
      printf("\nSequence file (MASE) ?  ");
      gets(nomfinseq);
      in=fopen(nomfinseq, "r");
      if(in) break;
      printf("Cannot find file : %s\n", nomfinseq);
    }
  }
  else{
    in=fopen(argv[1], "r");
    if(!in){
      printf("Cannot find sequence file : %s\n", argv[1]);
      exit(EXIT_FAILURE);
    }
    sprintf(nomfinseq, "%s", argv[1]);
    muet=1;
  }
  fclose(in);
  if(muet)
    nbseq=readmasemuet(nomfinseq, seq, seqname, comments, MAXDATASET*MAXNSP);
  else
    nbseq=readmaseseqs(nomfinseq, seq, seqname, comments, MAXDATASET*MAXNSP);
  nbseqvrai=nbseq;   
  for(i=1;i<=nbseq/2;i++){ 
    if(strncmp(comments[i], ";;", 2)==0){
      nbseqvrai=i;
      break;
    }
  }
  if(nbseq%nbseqvrai!=0){
    printf("Bad sequence file %s\n", nomfinseq);
    exit(EXIT_FAILURE);
  }
  nbdataset=nbseq/nbseqvrai;

  for(i=0;i<nbdataset;i++)
    refresh(seq+i*nbseqvrai, nbseqvrai, 0);


   /** input trees **/

  maxtreesize=50*nbseqvrai*MAXNTREE;
  alltrees=(char*)check_alloc(maxtreesize+1, sizeof(char));
  c_tree=(char**)check_alloc(MAXNTREE, sizeof(char*));

  if(argc<3){
    while(1){
      printf("\nTree file ?  ");
      gets(nomfintree);
      treefile=fopen(nomfintree, "r");
      if(treefile) break;
      printf("Cannot find file : %s\n", nomfintree);
    }
  }
  else{
    treefile=fopen(argv[2], "r");
    if(!treefile){
      printf("Cannot find tree file : %s\n", argv[2]);
      exit(EXIT_FAILURE);
    }
  }

  i=0;
  while(i<maxtreesize && (alltrees[i]=getc(treefile))!=EOF) i++;
  alltrees[i]='\0';

  fclose(treefile);

  i=0;
  prov=alltrees;
  while(*prov) {if(*prov==';') i++; prov++;}
  free(c_tree);
  c_tree=(char**)check_alloc(i, sizeof(char*));

  i=0; prov=alltrees;
  while(*prov){
    if(*prov=='[') i++;
    if(*prov==']') i--;
    if(i!=0 && i!=1){
      printf("Unmatched brackets [] in tree file\n");
      exit(EXIT_FAILURE);
    }
    if(i==1 && (*prov==';' || *prov=='(')) *prov='.';
    prov++;
  }


  nbtree=0;
  prov=alltrees;
  while(1){
    prov2=strtok(prov, ";");
    if(prov) prov=NULL;
    if(!prov2) break;
    while(*prov2 && *prov2!='(') prov2++;
    if(*prov2==0) continue;
    c_tree[nbtree]=prov2;
    nbtree++;
  }


  for(i=0;i<nbtree;i++){
    prov=c_tree[i];
    while(*prov) prov++;
    *prov=';'; *(prov+1)=0;
  }


   /** input options **/


  if(argc<4){
    while(1){
      printf("\nOption file ?  ");
      gets(nomfopt);
      optfile=fopen(nomfopt, "r");
      if(optfile) break;
      printf("Cannot find file : %s\n", nomfopt);
    }
  }
  else{
    optfile=fopen(argv[3], "r");
    if(!optfile){
      printf("Cannot find option file : %s\n", argv[3]);
      exit(EXIT_FAILURE);
    }
  }
  getoptions(&opt, optfile);
  print1=opt->print->PRINT1;


  allcouples=opt->ALLCOUPLES;

  if(opt->print->EVAL_OUT){
    dutheil=fopen("detailed_out", "w");
    fprintf(dutheil, "numero\tarbre\tlnL\tTs/Tv\tGCanc\talpha\tcovar\tpi\n");
    fclose(dutheil);
  }
  

  if(print1) 
    printf("\n%d sequence data sets and %d trees found : ", nbdataset, nbtree);
  if(!allcouples){
    nbdataset=mini(nbdataset, nbtree);
    nbtree=nbdataset;
  }
  if(print1)
    printf("%d evaluations processed\n\n", allcouples?nbdataset*nbtree:nbdataset);

   

  maxl=(double**)check_alloc(nbdataset, sizeof(double*));
  for(i=0;i<nbdataset;i++)
    maxl[i]=(double*)check_alloc(nbtree, sizeof(double));

  ctree1=(char*)check_alloc(50*nbseq, sizeof(char));
  ctree2=(char*)check_alloc(50*nbseq, sizeof(char));
  outfile1=fopen("treefile.eqgc", "w");
  outfile2=fopen("treefile.ndgc", "w");
  if(outfile1==NULL || outfile2==NULL){ printf("Cannot write tree file\n"); exit(EXIT_FAILURE); }

  for(i=0;i<nbdataset;i++){
    for(j=0;j<nbtree;j++){
      if(!allcouples) j=i;
      if(print1) printf("\ndata set %d , tree %d\n", i+1, j+1);
      maxl[i][j]=maxlike(nbseqvrai, seq+i*nbseqvrai, seqname+i*nbseqvrai, c_tree[j], opt, ctree1, ctree2);

      if(nbdataset*nbtree>1){
        fprintf(outfile1, "[data set %d, tree %d]\n%s\n\n", i+1, j+1, ctree1);
        fprintf(outfile2, "[data set %d, tree %d]\n%s\n\n", i+1, j+1, ctree2);
      }
      else{
        fprintf(outfile1, "%s\n", ctree1);
        fprintf(outfile2, "%s\n", ctree2);
      }
      if(!allcouples) break;
    }
  }
  
  if(print1){
    if(nbdataset*nbtree==1){
      printf("Tree is written into files : treefile.eqgc (equilibrium G+C content)\n");
      printf("                             treefile.ndgc (G+C content at each node)\n\n");
    }
    else{
      printf("Trees are written into files : treefile.eqgc (equilibrium G+C content)\n");
      printf("                               treefile.ndgc (G+C content at each node)\n\n");
    }
  }

}





