#include "nhmlg.h"

#define MIN_NBLISTMAX 1000
#define MAX_NBLISTMAX 10000


tree* s_tree;
compute_option compute_o;
noeud root;

extern int gamma_nbcl;

/*
char** debug_nom;
FILE* debugf;
int detailed_debug;

void debug_printree(int** ttree, int nb, int nbbi);
*/

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
void stoc(noeud* arbre_s, int racine, int notu, char* ctree, int nodegc);
void refresh(char** seq, int nbseq, int option);
double gg95(char* seq1, char* seq2, double alpha);
double compute(tree* s_tree, init_option init_o, compute_option compute_opt, converge_option converge_o, print_option print_o);
noeud create_node(noeud v1, noeud v2, noeud v3, double l1, double l2, double l3, double b1, double b2, double b3, char* nom);
int ctot(char *input, int *arbre[], double *lgbi, double *lgbp, double* bootstrap, char *nom[], char *racine, int* nbbi);
char* readtree(FILE* treef, int nbsp);
void getoptions(options* opt_arg, FILE* in_opt);





/* ctree_noblbs */
/* Copy c_tree in into c_tree out without branch length. */

void ctree_noblbs(char* in, char* out, int lgin){
  int i=0, ii=0, newi;
  char c;

  while(out[i]){out[i]=0; i++;}
  i=0;

  while(i<lgin){
    c=in[i];
    if (c!=':' && c!=' ') {out[ii]=c; ii++;}
    if (c==';') break;
    newi=i;
    if (c==')' || c==':' || c==' ') newi=i+1+strspn(in+i+1,"0123456789.-eE \n\t");
    if(newi!=i) i=newi; else i++;
  }
  out[ii]='\0';

}



noeud* ttos_unres_rooted(int** ttree, int nbotu, int nbbi, char** nom, double* lgbi_init, double* lgbp, double fracroot1, int l1, char** list1, int l2, char** list2, double** lgbi_pt){

  noeud *tabnd, nd1, nd2, nd3;
  int i, j, k, *kill_bi, min, rang_min, ** sorted, *nbfeuille, *fils, nbfils, cpt, *remain;
  double arg6, *lgbi;
  branche* br_list;

int mo1=-1;


  tabnd=(noeud*)check_alloc(2*nbotu, sizeof(noeud));
  kill_bi=(int*)check_alloc(nbbi, sizeof(int));
  remain=(int*)check_alloc(nbbi, sizeof(int));
  if(lgbi_init)
    lgbi=(double*)check_alloc(nbbi, sizeof(double));
  else lgbi=NULL;
  sorted=(int**)check_alloc(nbotu, sizeof(int*));
  for(i=0;i<nbotu;i++)
    sorted[i]=(int*)check_alloc(nbbi, sizeof(int));
  nbfeuille=(int*)check_alloc(nbbi, sizeof(int));
  fils=(int*)check_alloc(2*nbotu, sizeof(int));
  br_list=(branche*)check_alloc(2*nbotu-3, sizeof(branche));
  for(i=0;i<2*nbotu-3;i++){
    br_list[i]=(branche)check_alloc(1, sizeof(struct branche));
    br_list[i]->bout1=NULL;
  }

 
 		/* SORT BRANCHES , REFRESH TREE */

  for(i=0;i<nbbi;i++){
    nbfeuille[i]=0;
    for(j=0;j<nbotu;j++)
      if(ttree[j][i]==1) nbfeuille[i]++;
    if(nbfeuille[i]>nbotu/2){
      nbfeuille[i]=nbotu-nbfeuille[i];
      for(j=0;j<nbotu;j++) ttree[j][i]=1-ttree[j][i];
    }
  }

  for(i=0;i<nbbi;i++) kill_bi[i]=0;

  for(i=0;i<nbbi;i++){
    min=nbotu+1; rang_min=-1;
    for(j=0;j<nbbi;j++){
      if(nbfeuille[j]<min && kill_bi[j]==0){
	min=nbfeuille[j];
	rang_min=j;
      }
    }
    for(j=0;j<nbotu;j++) sorted[j][i]=ttree[j][rang_min];
    if(lgbi_init){ lgbi[i]=lgbi_init[rang_min]; remain[i]=rang_min; }
    kill_bi[rang_min]=1;
  }

   
  
 		/* CREATE NBOTU TERMINAL NODES */

  for(i=0;i<nbotu;i++){
    if(lgbp==NULL) arg6=-1.0; else arg6=lgbp[i];
    tabnd[i]=create_node(NULL, NULL, NULL, 0., 0., arg6, 0., 0., 1., nom[i]);
  }


  		/* CREATE NBBI ALIVE INTERNAL NODES + A FEW NONALIVE ONES*/

  cpt=0;
  for(i=0;i<nbbi;i++){

    		/* find and record sons */
    nbfils=0;
    for(j=0;j<nbbi;j++) kill_bi[j]=0;
    for(j=0;j<nbotu;j++){
      if(sorted[j][i]==1){
	for(k=i-1;k>=0;k--){
	  if(sorted[j][k]==1){
   	    if(kill_bi[k]) break;
/* internal son node (possible index : nbotu, nbotu+1, ..., nbotu+nbbi-1) */
	    nbfils++;
	    fils[nbfils-1]=nbotu+k;
	    kill_bi[k]=1;
	    break;
	  }
	}
	if(k==-1){
/* terminal son node (possible index : 0, 1, ..., nbotu-1) */
	  nbfils++;
	  fils[nbfils-1]=j;
	}
      }
    }
    if(nbfils<2) {
      printf("Bad tree:\n");
      printf("%d sons found for branch %d in sorted tree\n", nbfils, i);
      /* debug_printree(sorted, nbotu, nbbi);*/
      exit(0);
    }

		/* create node(s) */

    if(lgbi==NULL) arg6=-1.0; else arg6=lgbi[i];

    if(nbfils==2){
      nd1=tabnd[fils[0]];
      nd2=tabnd[fils[1]];
      tabnd[nbotu+i]=create_node(nd1, nd2, NULL, nd1->l3, nd2->l3, arg6, 1., 1., 1., NULL);
      nd1->v3=nd2->v3=tabnd[nbotu+i]; 
      if(lgbi_pt) lgbi_pt[remain[i]]=&(tabnd[nbotu+i]->l3);
    }
    else{
      nd1=tabnd[fils[0]]; nd2=tabnd[fils[1]];
      tabnd[nbotu+nbbi+cpt]=create_node(nd1, nd2, NULL, nd1->l3, nd2->l3, -1., 1., 1., 0., NULL);
      nd1->v3=nd2->v3=tabnd[nbotu+nbbi+cpt]; 
      cpt++;
      for(j=0;j<nbfils-3;j++){
	nd1=tabnd[fils[j+2]]; nd2=tabnd[nbotu+nbbi+cpt-1];
	tabnd[nbotu+nbbi+cpt]=create_node(nd1, nd2, NULL, nd1->l3, nd2->l3, -1., 1., 0., 0., NULL);
	nd1->v3=nd2->v3=tabnd[nbotu+nbbi+cpt]; 
	cpt++;
      }
      nd1=tabnd[fils[nbfils-1]]; nd2=tabnd[nbotu+nbbi+cpt-1];
      tabnd[nbotu+i]=create_node(nd1, nd2, NULL, nd1->l3, nd2->l3, arg6, 1., 0., 1., NULL);
      nd1->v3=nd2->v3=tabnd[nbotu+i]; 
    }


  }



		/* JOIN HANGING NODES */

  nbfils=0;
  for(i=0;i<nbotu+nbbi+cpt;i++){
    if(tabnd[i]->v3==NULL){
      fils[nbfils]=i;
      nbfils++;
    }
  }

 
  if(nbfils<3){
    printf("Bad tree 2\n");
    exit(0);
  }

  if(nbfils==3){
    nd1=tabnd[fils[0]];
    nd2=tabnd[fils[1]];
    nd3=tabnd[fils[2]];
    tabnd[nbotu+nbbi+cpt]=create_node(nd1, nd2, nd3, nd1->l3, nd2->l3, nd3->l3, 1., 1., 1., NULL);
    nd1->v3=nd2->v3=nd3->v3=tabnd[nbotu+nbbi+cpt];
  }

  if(nbfils>3){
    nd1=tabnd[fils[0]]; nd2=tabnd[fils[1]];
    tabnd[nbotu+nbbi+cpt]=create_node(nd1, nd2, NULL, nd1->l3, nd2->l3, -1., 1., 1., 0., NULL);
    nd1->v3=nd2->v3=tabnd[nbotu+nbbi+cpt];
    cpt++;
    for(i=0;i<nbfils-4;i++){
      nd1=tabnd[fils[i+2]]; nd2=tabnd[nbotu+nbbi+cpt-1];
      tabnd[nbotu+nbbi+cpt]=create_node(nd1, nd2, NULL, nd1->l3, nd2->l3, -1., 1., 0., 0., NULL);
      nd1->v3=nd2->v3=tabnd[nbotu+nbbi+cpt];
      cpt++;
    }
    nd1=tabnd[fils[nbfils-2]]; nd2=tabnd[fils[nbfils-1]]; nd3=tabnd[nbotu+nbbi+cpt-1];
    tabnd[nbotu+nbbi+cpt]=create_node(nd1, nd2, nd3, nd1->l3, nd2->l3, nd3->l3, 1., 1., 0., NULL);
    nd1->v3=nd2->v3=nd3->v3=tabnd[nbotu+nbbi+cpt];
  }
 		/* ROOT */



  root=create_node(NULL, NULL, NULL, -1., -1., -1., 1., 1., 0., "ROOT");

  makelistbr_unrooted(NULL, tabnd[0], br_list, NULL, 2*nbotu-3);


  for(i=0;i<2*nbotu-3;i++){
    if(isbranch(br_list[i], l1, list1, l2, list2)==1){
      root_s(br_list[i], fracroot1);
      break;
    }
    if(isbranch(br_list[i], l1, list1, l2, list2)==2){
      root_s(br_list[i], 1.-fracroot1);
      break;
    }
  }
  if(i==2*nbotu-3) printf("Erreur racinage\n");

  organize_tree(NULL, root);


  free(kill_bi);
  if(lgbi_init) free(lgbi);
  for(i=0;i<nbotu;i++) free(sorted[i]);
  free(sorted);
  free(nbfeuille); free(fils);
  for(i=0;i<2*nbotu-3;i++) free(br_list[i]);
  free(br_list);
  free(remain);


  return(tabnd);
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

/*
void debug_printclade(int** ttree, int nb, int nbbi){

  int i, j;

  fprintf(debugf, "       ");
  for(i=0;i<nbbi;i++) fprintf(debugf, "%d", i%10);
  fprintf(debugf, "\n");

  for(i=0;i<nb;i++){
    fprintf(debugf, "%d: ", i%10);
    for(j=0;j<nbbi;j++) fprintf(debugf, "%d", ttree[i][j]);
    fprintf(debugf, "\n");
  }
  fprintf(debugf, "\n");
}



void debug_printree(int** ttree, int nb, int nbbi){

  int i, j;

  fprintf(debugf, "       ");
  for(i=0;i<nbbi;i++) fprintf(debugf, "%d", i%10);
  fprintf(debugf, "\n");

  for(i=0;i<nb;i++){
    fprintf(debugf, "%.5s: ", debug_nom[i]);
    for(j=0;j<nbbi;j++) fprintf(debugf, "%d", ttree[i][j]);
    fprintf(debugf, "\n");
  }
  fprintf(debugf, "\n");
}
*/


double maxlike(int nbseq, char** seq, char** seqname, int** ttree, double* lgbi, double* lgbp, int nbbi, int l1, char** list1, int l2, char** list2, options opt, char* ctree, char* ctree1, char* ctree2){

  int i, j, k, l, ii, lgseq, lgsite, *drapeau, *drapeau2, cont=0, nbparam, nbvrai;
  char **orderedseq, racine, *site;
  double *weight, lkh, maxlkh=-1.e100, prov, **lgbi_pt, *freq;

  gamma_nbcl=opt->init->GAMMA_NCL;

  /*  INITIALIZATION */

  s_tree=(tree*)check_alloc(1, sizeof(tree));
  lgseq=(int)strlen(seq[0]);
  

  s_tree->alivebr=(int*)check_alloc(2*nbseq, sizeof(int));
  if(lgbi) 
    lgbi_pt=(double**)check_alloc(nbseq-2, sizeof(double*));
  else
    lgbi_pt=NULL;
  weight=(double*)check_alloc(lgseq, sizeof(double));
/*
debug_printree(ttree, nbseq, nbbi);
*/

  s_tree->node=ttos_unres_rooted(ttree, nbseq, nbbi, seqname, lgbi, lgbp, 0.5, l1, list1, l2, list2, lgbi_pt);

  s_tree->nbseq=nbseq;

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


  for(i=0;i<2*nbseq-2;i++)
    alloc_node(s_tree->node[i], lgseq, nbparam, nbseq, (i<nbseq), 0, gamma_nbcl);
  alloc_node(root, lgseq, nbparam, nbseq, 0, 0, gamma_nbcl);


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

  for(i=nbseq;i<2*nbseq-2;i++) sprintf(s_tree->node[i]->nom, "int%d", i-nbseq+1);

  drapeau=(int*)check_alloc(nbseq, sizeof(int));
  drapeau2=(int*)check_alloc(nbseq, sizeof(double));
  orderedseq=(char**)check_alloc(nbseq, sizeof(char*));
  for(i=0;i<nbseq;i++) drapeau2[i]=0;
  

  for(i=0;i<nbseq;i++){
    for(j=0;j<nbseq;j++){
      if(!drapeau[j] && strncmp(s_tree->node[i]->nom, seqname[j], MAXLNAME)==0){
        drapeau[j]=drapeau2[i]=1;
	orderedseq[i]=seq[j];
        break;
      }
    }
    if(drapeau2[i]==0){
      printf("Taxon %s not found in sequence file.\n", s_tree->node[i]->nom);
      exit(0);
    }
  }
  s_tree->nbseq=nbseq;

  for(i=0;i<lgseq;i++) weight[i]=0.;
  site=(char*)check_alloc(nbseq+1, sizeof(char)); 
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
    if(cont) continue; 
    for(j=0;j<nbseq;j++)
      s_tree->node[j]->seq[lgsite]=site[j];
    weight[lgsite]=1.;
    lgsite++;
  }


  s_tree->lgseq=lgsite;
  s_tree->weight=weight;
  s_tree->weightsum=lgseq;


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
    opt->init->INIT_TITV/=nbvrai;
    free(freq);
  }

  if(opt->init->INIT_ROOT<=0.){
    if(lgbp[0]>-0.5)
      opt->init->INIT_ROOT=root->l1/(root->l1+root->l2);
    else
      opt->init->INIT_ROOT=0.5;
  }

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


  organize_tree(NULL, root);
  organize_listbr(s_tree);
  setpattern(s_tree);
  for(ii=0;ii<2*nbseq-3;ii++)
    setunderlyingbr(root, s_tree->listbr[ii], ii);
  setunderlyinggc(root, s_tree->nbseq);
  set_numerobr(root, s_tree->listbr, 2*nbseq-3);



  lkh=compute(s_tree, opt->init, opt->compute, opt->converge, opt->print);
  if(lkh>maxlkh) maxlkh=lkh;

  if(opt->print->PRINT0)
    printf("ln(L)=%f\n", lkh);


  if(ctree)
    stoc(s_tree->node, 1, nbseq, ctree, -1);

  if(ctree1)
    stoc(s_tree->node, 1, nbseq, ctree1, 0);

  if(ctree2){
    set_node_deduced_comp(root);
    stoc(s_tree->node, 1, nbseq, ctree2, 1);
  }



  if(lgbi){
    for(i=0;i<nbseq-3;i++){
      lgbi[i]=*(lgbi_pt[i]);
    }
  }



  /* FREEING */

  free(drapeau); free(drapeau2); free(site); free(orderedseq);
  if(lgbi) free(lgbi_pt);

  free_tree(s_tree, gamma_nbcl);
  free(s_tree);
  free_node(root, lgseq, nbparam, 0, gamma_nbcl);


  return maxlkh;

}



void setnewbranch(int** tree, int nbotu, int numbi, int* listotu, int* listotu_init, int listinitsize, int j, int k){
  int l;

  for(l=0;l<nbotu;l++) tree[l][numbi]=0;
  for(l=0;l<listinitsize;l++)
    if(abs(listotu[l])==listotu[j] || abs(listotu[l])==listotu[k])
      tree[listotu_init[l]][numbi]=1;
}



/* Writes in ranks the ranks of tab values (0 -> nb-1) */
/* Destroys tab */

void sort(int* tab, int* ranks, int nb){
  int i, j, max, rgmax;

  
  for(i=0;i<nb;i++){
    max=rgmax=-1;
    for(j=0;j<nb;j++){
      if(tab[j]>max){
        max=tab[j];
	rgmax=j;
      }
    }
    if(rgmax==-1) break;
    ranks[rgmax]=nb-1-i;
    tab[rgmax]=-2;
  }

}




void nesteddist(int** clade, int nb, int nbcl, int* clsz, int otu, int** dist){
  int i, j, k, *cllist, *clldist, *ranks, nbnest=0;

  cllist=(int*)check_alloc(nbcl, sizeof(int));
  clldist=(int*)check_alloc(nbcl, sizeof(int));
  ranks=(int*)check_alloc(nbcl, sizeof(int));


  for(i=0;i<nbcl;i++){
    dist[i][i]=0;
    if(clade[i][otu]){
      cllist[nbnest]=i;
      clldist[nbnest]=clsz[i];
      nbnest++;
    }
  }

  sort(clldist, ranks, nbnest);

  for(j=0;j<nbnest;j++){
    for(k=j+1;k<nbnest;k++){
      dist[cllist[j]][cllist[k]]=ranks[k]-ranks[j];
      if(dist[cllist[j]][cllist[k]]>0)
	dist[cllist[j]][cllist[k]]--;
      else
	dist[cllist[j]][cllist[k]]++;
      dist[cllist[k]][cllist[j]]=-dist[cllist[j]][cllist[k]];
    }
  }

/* For nested clades, dist[i][j]=-dist[j][i]. dist[i][j] is positive if     */
/* clade i is a part of clade j, negative if clade j is a part of clade i,  */
/* so that d[i][j]<0 => clade i cannot move toward clade j.		    */
/* Non-nested clades have positive distances both ways. */

  free(cllist); free(clldist); free(ranks);
}



/* Returns the rank of the last common ancestral clade of non-nested clades cl1 and cl2 */

int lca(int** clade, int nbcl, int nbotu, int cl1, int cl2){

  int i, j, *lcacl, noparent, sz, min, rgmin;

  lcacl=(int*)check_alloc(nbotu, sizeof(int));
  
  for(i=0;i<nbotu;i++){
    if(clade[cl1][i]==1 || clade[cl2][i]==1) lcacl[i]=1;
    else  lcacl[i]=0;
  }
  
  min=nbotu+1; rgmin=-1;
  for(i=0;i<nbcl;i++){
    sz=0; noparent=0;
    for(j=0;j<nbotu;j++){
      if(clade[i][j]==0 && lcacl[j]==1) {noparent=1; break;}
      sz+=clade[i][j];
    }
    if(noparent) continue;
    if(sz<min){
      min=sz;
      rgmin=i;
    }
  }

  if(rgmin<0){
    printf("erreur lca\n");
    exit(0);
  }

  free(lcacl);

  return rgmin;
}


void  setclades(int** tree, int nbotu, char** nom, int l1, char** list1, int l2, char** list2, int** dclade, int** gclade, int**ddist, int** gdist, int nbdclade, int nbgclade, char** dcladename, char** gcladename){

  int i, j, k, l, d0, d1, g0, g1, *prov, **prov2, ** dprov2, nbd=0, nbg=0, sum, val;
  int *cladesz;
  char c1, c2, **provc, **clname, *ch, *droitgauche;


  droitgauche=(char*)check_alloc(nbotu, sizeof(char));

  if(l1>l2) {k=l2; provc=list2; c1='g'; c2='d';}
  else {k=l1; provc=list1; c1='d'; c2='g';}

  for(i=0;i<nbotu;i++){
    for(j=0;j<k;j++){
      if(samename(nom[i], provc[j])){
	droitgauche[i]=c1;
	break;
      }
    }
    if(j==k) droitgauche[i]=c2;
  }

		/* peripheric clades */

  for(i=0;i<nbotu;i++){
    if(droitgauche[i]=='d') {prov2=dclade; prov=&nbd; clname=dcladename;} 
    else { prov2=gclade; prov=&nbg; clname=gcladename;}
    for(j=0;j<nbotu;j++)
      prov2[*prov][j]=0;
    prov2[*prov][i]=1;
    sprintf(clname[*prov], "(%s)", nom[i]);
    *prov=*prov+1;
  }


		/* internal clades */

  for(i=nbd;i<nbdclade;i++){
    sprintf(dcladename[i], "(");
  }
  for(i=nbg;i<nbgclade;i++){
    sprintf(gcladename[i], "(");
  }

  for(i=0;i<nbotu-3;i++){
    d0=d1=g0=g1=0;
    for(j=0;j<nbotu;j++){
      if(tree[j][i]==0 && droitgauche[j]=='d') d0=1;
      if(tree[j][i]==1 && droitgauche[j]=='d') d1=1;
      if(tree[j][i]==0 && droitgauche[j]=='g') g0=1;
      if(tree[j][i]==1 && droitgauche[j]=='g') g1=1;
    }
    sum=d0+d1+g0+g1;
    if(sum!=2 && sum!=3) {printf("probleme faire clades\n"); exit(0);}

    if(sum==2){ /* root branch -> 2 full clades */
      for(j=0;j<nbotu;j++)
	if(droitgauche[j]=='d'){
	  dclade[nbd][j]=1;
	  strcat(dcladename[nbd], nom[j]);
	  strcat(dcladename[nbd], ",");
	}
      for(j=0;j<nbotu;j++) 
	if(droitgauche[j]=='g'){ 
	  gclade[nbg][j]=1;
	  strcat(gcladename[nbg], nom[j]);
	  strcat(gcladename[nbg], ",");
	}
      ch=dcladename[nbd];
      while(*ch) ch++; while(*ch!=',' && *ch!=')') ch--; *ch=')';
      ch=gcladename[nbg];
      while(*ch) ch++; while(*ch!=',' && *ch!=')') ch--; *ch=')';
      nbd++; nbg++;
    }
    else if(d0 && d1){  /* right clade */
      for(j=0;j<nbotu;j++){
	if(g0) val=tree[j][i];
	else val=1-tree[j][i];
	dclade[nbd][j]=val;
	if(val){
	  strcat(dcladename[nbd], nom[j]);
	  strcat(dcladename[nbd], ",");
	}
      }
      ch=dcladename[nbd];
      while(*ch) ch++; while(*ch!=',' && *ch!=')') ch--; *ch=')';
      nbd++;
    }
    else if(g0 && g1){  /* left clade */
      for(j=0;j<nbotu;j++){
	if(d0) val=tree[j][i];
	else val=1-tree[j][i];
	gclade[nbg][j]=val;
	if(val){
	  strcat(gcladename[nbg], nom[j]);
	  strcat(gcladename[nbg], ",");
	}
      }
      ch=gcladename[nbg];
      while(*ch) ch++; while(*ch!=',' && *ch!=')') ch--; *ch=')';
      nbg++;
    }
    else{ printf("erreur setting clades\n"); exit(0); }
  }

		/* last clade */

  if(nbd<nbdclade){
    for(i=0;i<nbotu;i++){
      if(droitgauche[i]=='d'){
	dclade[nbd][i]=1;
	strcat(dcladename[nbd], nom[i]);
	strcat(dcladename[nbd], ",");
      }
      else
	dclade[nbd][i]=0;
    }
    ch=dcladename[nbd];
    while(*ch) ch++; while(*ch!=',' && *ch!=')') ch--; *ch=')';
    nbd++;
  }

  if(nbg<nbgclade){
    for(i=0;i<nbotu;i++){
      if(droitgauche[i]=='g'){
	gclade[nbg][i]=1;
	strcat(gcladename[nbg], nom[i]);
	strcat(gcladename[nbg], ",");
      }
      else
	gclade[nbg][i]=0;
    }
    ch=gcladename[nbg];
    while(*ch) ch++; while(*ch!=',' && *ch!=')') ch--; *ch=')';
    nbg++;
  }

  if(nbd!=nbdclade || nbg!=nbgclade) {
    printf("wrong number of clades\n");
    exit(0);
  }
  		/* distances */

  if(nbd>nbg) k=nbd; else k=nbg;
  cladesz=(int*)check_alloc(k, sizeof(int));


  for(k=0;k<2;k++){
    if(k==0){ 		/* right */
      prov=&nbd;
      prov2=dclade;
      for(i=0;i<*prov;i++)
        for(j=0;j<*prov;j++)
	  ddist[i][j]=-2*nbotu;
      dprov2=ddist;
      c1='d';
    }
    else{		/* left  */
      prov=&nbg;
      prov2=gclade;
      for(i=0;i<*prov;i++)
        for(j=0;j<*prov;j++)
	  gdist[i][j]=-2*nbotu;
      dprov2=gdist;
      c1='g';
    }

  /* sizes of clades */
    for(i=0;i<*prov;i++){
      cladesz[i]=0;
      for(j=0;j<nbotu;j++) 
	cladesz[i]+=prov2[i][j]; 
    }

  /* distances between nested clades */
    for(i=0;i<*prov;i++){
      if(droitgauche[i]==c1)
        nesteddist(prov2, nbotu, *prov, cladesz, i, dprov2);  
    }

  /* distances between non-nested clades */
    for(i=0;i<*prov;i++){
      for(j=i+1;j<*prov;j++){
	if(dprov2[i][j]>-2*nbotu) continue;
	l=lca(prov2, *prov, nbotu, i, j);
	dprov2[i][j]=dprov2[j][i]=abs(dprov2[i][l])+abs(dprov2[j][l]);
      }
    }
  }

  free(droitgauche);
  free(cladesz);
}



void moveclade(int** clade, int** newclade, int nbclade, int nbotu, int** dist, int cl, int tocl){

  int i, j, k, cl_a_virer;
  int one_in, one_out, zero_in;

	/* copy clade into newclade */

  for(i=0;i<nbclade;i++){
    for(j=0;j<nbotu;j++)
      newclade[i][j]=clade[i][j];
  }
/*
if(detailed_debug){
  fprintf(debugf, "clades:\n");
  debug_printclade(clade, nbclade, nbotu);
}
*/

	/* find and modify clades affected by the departure of cl */

  for(i=0;i<nbclade;i++){
    one_in=one_out=0;
    for(j=0;j<nbotu;j++){
      if(clade[i][j] && clade[cl][j]) one_in=1;
      if(clade[i][j] && !clade[cl][j]) one_out=1;
      if(one_in && one_out) break;
    }
    if(one_in && one_out){
      for(j=0;j<nbotu;j++) if(clade[cl][j]) newclade[i][j]=0;
    }
  }

/*
if(detailed_debug){
  fprintf(debugf, "clades after departure:\n");
  debug_printclade(newclade, nbclade, nbotu);
}
*/

	/* determine and store the redundant clade */

  cl_a_virer=-1;
  for(i=0;i<nbclade;i++){
    for(j=i+1;j<nbclade;j++){
      for(k=0;k<nbotu;k++){ 
	if(newclade[i][k]!=newclade[j][k]){
	  break;
	}
      }
      if(k==nbotu){
	cl_a_virer=j;
	break;
      }
    }
    if(cl_a_virer>-1) break;
  }
/*
fprintf(debugf, "cl=%d, tocl=%d, cl_a_virer=%d\n", cl, tocl, cl_a_virer);
*/
  if(cl_a_virer<0){
    printf("erreur red\n");
    exit(0);
  }

	/* find and modify clades affected by the arrival of cl */


  for(i=0;i<nbclade;i++){
    if(i==cl || i==tocl) continue;
    zero_in=0;
    for(j=0;j<nbotu;j++){
      if(clade[i][j]==0 && clade[tocl][j]==1){
	zero_in=1;
	break;
      }
    }
    if(!zero_in){
      for(j=0;j<nbotu;j++)
	if(clade[cl][j]) newclade[i][j]=1;
    }
  }

/*
if(detailed_debug){
  fprintf(debugf, "clades after arrival:\n");
  debug_printclade(newclade, nbclade, nbotu);
}
*/


	/* replace the redundant clade by new (cl, tocl) clade */

  for(j=0;j<nbotu;j++){
    if(clade[cl][j] || clade[tocl][j])
      newclade[cl_a_virer][j]=1;
    else
      newclade[cl_a_virer][j]=0;
  }

/*
if(detailed_debug){
  fprintf(debugf, "clades fin:\n");
  debug_printclade(newclade, nbclade, nbotu);
}
*/

}



int samebranch(int* cl1, int* cl2, int nbotu){

  int i, eq, opp;

  eq=0; opp=0;
  for(i=0;i<nbotu;i++){
    if(cl1[i]==cl2[i]) eq=1;
    else opp=1;
    if(eq && opp) return 0;
  }
  return 1;
}



void settree(int** dclade, int nbdclade, int** gclade, int nbgclade, int nbotu, int** newtree){

  int i, j, k, a_virer_g=-1, sz, nbbi=0;

  for(i=0;i<nbdclade;i++){
    for(j=0;j<nbgclade;j++){
      if(samebranch(dclade[i], gclade[j], nbotu)){
	a_virer_g=j;
	break;
      }
    }
    if(a_virer_g>-1) break;
  }

  for(i=0;i<nbdclade;i++){
    sz=0;
    for(j=0;j<nbotu;j++){
      sz+=dclade[i][j];
      if(sz>=2) break;
    }
    if(j!=nbotu){
      for(k=0;k<nbotu;k++)
	newtree[k][nbbi]=dclade[i][k];
      nbbi++;
    }
  }

  for(i=0;i<nbgclade;i++){
    if(i==a_virer_g) continue;
    sz=0;
    for(j=0;j<nbotu;j++){
      sz+=gclade[i][j];
      if(sz>=2) break;
    }
    if(j!=nbotu){
      for(k=0;k<nbotu;k++)
	newtree[k][nbbi]=gclade[i][k];
      nbbi++;
    }
  }

}



int deja_evalue(int** newtree, int nb, int** list_tree, int nblist, int nblistmax){

  int i, ii, j, jj, k, debut, tot, bitdeb;
  int eq, opp, ret;



  debut=nblist%nblistmax;
  if(nblist<nblistmax) tot=nblist; else tot=nblistmax;
  for(i=0;i<tot;i++){
    ii=debut-i-1;
    if(ii<0) ii+=nblistmax;
    bitdeb=ii*(nb-3);
    for(j=0;j<nb-3;j++){
      for(jj=0;jj<nb-3;jj++){
        eq=opp=0;
        for(k=0;k<nb;k++){
	  ret=testbit(list_tree[k], bitdeb+jj+1);
	  if(ret && newtree[k][j]) eq=1;
	  else if(ret && !newtree[k][j]) opp=1;
	  else if(!ret && newtree[k][j]) opp=1;
	  else eq=1;
	  if(eq && opp) break; /* branches j (newtree) et jj (listed current tree) differentes */
        }
        if(k==nb) break; /* branches j et jj identiques : passer a j+1 */
      }
      if(jj==nb-3) break; /* branche j n'a pas d'identique : passer au listed tree suivant */
    }
    if(j==nb-3) return 1; /* toutes les branches j ont un identique : deja evalue */ 
  }
  return 0;
}



int long_branch_crossed(int** newtree, int** curtree, int nb, double* lgbi, double maxlcrossedbranch){

  int i, j, k, eq, opp, l10, l11, l20, l21;

  for(i=0;i<nb-3;i++){
    for(j=0;j<nb-3;j++){
      eq=opp=0;
      for(k=0;k<nb;k++){
	if(curtree[k][i] && newtree[k][j]) eq=1;
	else if(curtree[k][i] && !newtree[k][j]) opp=1;
	else if(!curtree[k][i] && newtree[k][j]) opp=1;
	else eq=1;
	if(eq && opp) break; 	
      }
      if(k==nb)
	break;   
    }
    if(j==nb-3){	/* internal branch i in tree curtree is absent in tree newtree */
      if(lgbi[i]>maxlcrossedbranch) 
 	return 1;
    }
  }

  return 0;

}


int solid_branch_crossed(int** newtree, int** curtree, int nb, int* solid){

  int i, j, k, eq, opp, l10, l11, l20, l21;

  for(i=0;i<nb-3;i++){
    for(j=0;j<nb-3;j++){
      eq=opp=0;
      for(k=0;k<nb;k++){
	if(curtree[k][i] == newtree[k][j]) eq=1;
	else opp=1;
	if(eq && opp) break; 	
      }
      if(k==nb)
	break;
    }
    if(j==nb-3){	/* internal branch i in tree curtree is absent in tree newtree */
      if(solid[i])
 	return 1;
    }
  }

  return 0;

}




void addtolist(int** newtree, int nb, int** list_tree, int nblist, int nblistmax){

  int i, j, debut;

  debut=nblist%nblistmax;
  for(i=0;i<nb;i++){
    for(j=0;j<nb-3;j++){
      if(newtree[i][j]) bit1(list_tree[i], debut*(nb-3)+j+1);
    }
  }

}


void copytree(int** from, int** to, int nb){

  int i, j, jfrom, jto;
  int eq, opp;
  int *keep_to, *keep_from;

  keep_to=check_alloc(nb, sizeof(int));
  keep_from=check_alloc(nb, sizeof(int));

/* seek common branches */

  for(jfrom=0;jfrom<nb-3;jfrom++){
    for(jto=0;jto<nb-3;jto++){
      if(keep_to[jto]==1) continue;
      eq=opp=0;
      for(i=0;i<nb;i++){
        if(from[i][jfrom]==to[i][jto]) eq=1;
        else opp=1;
        if(eq && opp) break;
      }
      if(i==nb){ /* identical branches */
        keep_to[jto]=1;
        keep_from[jfrom]=1;
        break;
      }
    }
  }


/* copy distinct branches */

  for(jfrom=0;jfrom<nb-3;jfrom++){
    if(keep_from[jfrom]==1) continue;
    for(jto=0;jto<nb-3;jto++){
      if(keep_to[jto]==1) continue;
      for(i=0;i<nb;i++)
        to[i][jto]=from[i][jfrom];
      keep_to[jto]=1;
      break;
    }
  }

  free(keep_to);
  free(keep_from);

}



void save_current_best(char* nom_f, char* ctree, double lkh, int accu){

  FILE* current_best;

  current_best=fopen(nom_f, "w");
  fprintf(current_best, "%s %f", ctree, lkh);
  if(accu) fprintf(current_best, "(ACCURATE)");
  fprintf(current_best, "\n");
  fclose(current_best);

}



void save_evaluated(char* nom_f, char* ctree, double lkh, int accu){
  FILE* evaluated;

  evaluated=fopen(nom_f, "a");
  fprintf(evaluated, "%s %f", ctree, lkh);
  if(accu) fprintf(evaluated, " (ACCURATE)");
  fprintf(evaluated, "\n");
  fclose(evaluated);

}



double shake(int nb, char** seq, char** seqname, char* ctree, options opt, char** eval_input, int nb_eval_input){

  int nb2, i, j, ii, jj, k, l, *ttree[MAXNSP], **curtree, **newtree, **evaltree, nbbi, nbdclade, nbgclade, print1, print2;
  int nochange, pres_grossiere=-1, l1, l2, restart_d, restart_g, oldmovedist;
  int **dclade, **gclade, **newdclade, **newgclade, **ddist, **gdist, movedist, maxmovedist;
  int **list_tree, lliste, *solid;
  long nblist, nblistmax;
  double *lgbp, *sortedlgbp, *lgbi, *bootvals, fracroot1, lkh, maxlkh, maxlcrossedbranch;
  char *nom[MAXNSP], *nom2[MAXNSP], **dcladename, **gcladename, **list1, **list2, racine, *ctreenew, *ctreenew_nobl, *treedeb, *ctree1, *ctree2;
  FILE* outfile1, *outfile2;
  print_option trueprint, noprint;


  print1=opt->print->PRINT1;
  print2=opt->print->PRINT2;
  noprint=check_alloc(1, sizeof(struct print_option));
  noprint->PRINT3=0;
  if(opt->print->PRINT2)
    noprint->PRINT1=noprint->PRINT2=1;
  else
    noprint->PRINT1=noprint->PRINT2=0;
  if(print1) noprint->PRINT0=1; else noprint->PRINT0=0;
  trueprint=opt->print;
  opt->print=noprint;
  maxlcrossedbranch=opt->SH_MAXLCROSSED;



		/* READ TREE STRING */

  lgbp=(double*)check_alloc(nb+1, sizeof(double));
  sortedlgbp=(double*)check_alloc(nb+1, sizeof(double));
  lgbi=(double*)check_alloc(nb+1, sizeof(double));
  for(i=0;i<nb+1;i++)
    lgbp[i]=lgbi[i]=-1.;

  solid=(int*)check_alloc(nb+1, sizeof(int));
  bootvals=(double*)check_alloc(nb+1, sizeof(double));
  if(nb>=MAXNSP) {printf("Too many sequences\n"); exit(0);}
  for(i=0;i<=nb;i++){
    nom[i]=(char*)check_alloc(MAXLNAME+1, sizeof(char));
    nom2[i]=(char*)check_alloc(MAXLNAME+1, sizeof(char));
    ttree[i]=(int*)check_alloc(nb, sizeof(int));
  }
  list1=check_alloc(nb, sizeof(char*));
  list2=check_alloc(nb, sizeof(char*));
  curtree=(int**)check_alloc(nb+1, sizeof(int*));
  newtree=(int**)check_alloc(nb+1, sizeof(int*));
  evaltree=(int**)check_alloc(nb+1, sizeof(int*));
  for(i=0;i<nb+1;i++) newtree[i]=(int*)check_alloc(nb-2, sizeof(int));
  ctreenew=(char*)check_alloc(2*MAXLNAME*nb, sizeof(char));
  ctree1=(char*)check_alloc(2*MAXLNAME*nb, sizeof(char));
  ctree2=(char*)check_alloc(2*MAXLNAME*nb, sizeof(char));
  ctreenew_nobl=(char*)check_alloc(2*MAXLNAME*nb, sizeof(char));


  nb2=ctot(ctree, ttree, lgbi, lgbp, bootvals, nom, &racine, &nbbi);

  if(nb2<nb){
    printf("More species in sequence file than in tree file\n");
    exit(0);
  }
  if(nb2>nb){
    printf("More species in tree file than in sequence file\n");
    exit(0);
  }


	/* PREPARE TREE : UNROOT, SET LEFT and RIGHT LISTS, SORT TAXA */

  if(racine=='r'){
    unroot(ttree, nb, lgbi, lgbp, bootvals, nom, list1, list2, &l1, &l2, &fracroot1);
  }
  else{
    printf("Tree must be rooted\n");
    exit(0);
  }
  nbbi--;
  if(nbbi!=nb-3){
    printf("Tree must be bifurcating\n");
    exit(0);
  }

  for(i=0;i<nb;i++){
    for(j=0;j<nb;j++){
      if(samename(seqname[i], nom[j])){
	curtree[i]=ttree[j];
	sortedlgbp[i]=lgbp[j];
	break;
      }
    }
  }

  for(i=0;i<nb-3;i++) if(bootvals[i]>opt->SH_MAXBOOTCROSSED) solid[i]=1;

  ctree_noblbs(ctree, ctreenew_nobl, strlen(ctree));

	/* EVALUATE INITIAL TREE */

  if(print1)
    printf("\nEvaluating initial tree : \n%s\n", ctreenew_nobl);
  if(print2)
    printf("\n");

  maxlkh=maxlike(nb, seq, seqname, curtree, lgbi, sortedlgbp, nbbi, l1, list1, l2, list2, opt, NULL, NULL, NULL);
  if(opt->SH_RESTART) save_current_best("current_best_tree", ctreenew_nobl, maxlkh, 1);
  if(opt->SH_RESTART) save_evaluated("evaluated_trees", ctreenew_nobl, maxlkh, 1);

	/* ALLOCATE SHAKE VARIABLES */

  nbdclade=2*l1-1;
  nbgclade=2*l2-1;
  dclade=(int**)check_alloc(nbdclade, sizeof(int*));
  for(i=0;i<nbdclade;i++)
    dclade[i]=(int*)check_alloc(nb, sizeof(int));
  gclade=(int**)check_alloc(2*nb, sizeof(int*));
  for(i=0;i<nbgclade;i++)
    gclade[i]=(int*)check_alloc(nb, sizeof(int));
  newdclade=(int**)check_alloc(nbdclade, sizeof(int*));
  for(i=0;i<nbdclade;i++)
    newdclade[i]=(int*)check_alloc(nb, sizeof(int));
  newgclade=(int**)check_alloc(nbgclade, sizeof(int*));
  for(i=0;i<nbgclade;i++)
    newgclade[i]=(int*)check_alloc(nb, sizeof(int));
  ddist=(int**)check_alloc(nbdclade, sizeof(int*));
  gdist=(int**)check_alloc(nbgclade, sizeof(int*));
  for(i=0;i<nbdclade;i++) ddist[i]=(int*)check_alloc(nbdclade, sizeof(int));
  for(i=0;i<nbgclade;i++) gdist[i]=(int*)check_alloc(nbgclade, sizeof(int));
  dcladename=(char**)check_alloc(nbdclade, sizeof(char*));
  gcladename=(char**)check_alloc(nbgclade, sizeof(char*));
  for(i=0;i<nbdclade;i++)
    dcladename[i]=(char*)check_alloc(nb*(MAXLNAME+3)+1, sizeof(char));
  for(i=0;i<nbgclade;i++)
    gcladename[i]=(char*)check_alloc(nb*(MAXLNAME+3)+1, sizeof(char));
  nblistmax=nb*nb;
  if(nblistmax<MIN_NBLISTMAX) nblistmax=MIN_NBLISTMAX;
  if(nblistmax<nb_eval_input) nblistmax=nb_eval_input;
  if(nblistmax>MAX_NBLISTMAX) nblistmax=MAX_NBLISTMAX;
  while(1){
    lliste=(nblistmax*(nb-3)+lmot-1)/lmot;
    list_tree=(int**)check_alloc(nb, sizeof(int*));
    for(i=0;i<nb;i++)
      list_tree[i]=(int*)calloc(lliste, sizeof(int));
    if(list_tree[nb-1]) break;
    nblistmax/=2;
    if(nblistmax==0){
      printf("Not enough memory\n");
      exit(0);
    }
  }


	/* SET LIST OF EVALUATED TREES */

  if(eval_input){
    nblist=0;
    for(k=0;k<nb_eval_input;k++){
      for(i=0;i<=nb;i++) ttree[i]=check_alloc(nb, sizeof(int));
      nb2=ctot(eval_input[k], ttree, NULL, NULL, NULL, nom2, &racine, NULL);
      if(racine=='r')
        unroot(ttree, nb, NULL, NULL, NULL, nom2, NULL, NULL, NULL, NULL, NULL);
      else{ printf("Evaluated trees must be rooted\n"); exit(0);}
      for(i=0;i<nb;i++){
        for(j=0;j<nb;j++){
          if(samename(seqname[i], nom2[j])){
	    evaltree[i]=ttree[j];
	    break;
          }
        }
      }
      if(deja_evalue(evaltree, nb, list_tree, nblist, nblistmax)) continue;
      addtolist(evaltree, nb, list_tree, nblist, nblistmax);
      nblist++;
    }
    printf("%d already evaluated topologies loaded\n", nblist);
  }
  else{
    addtolist(curtree, nb, list_tree, 0, nblistmax);
    nblist=1;
  }


	/* SHAKE */

  treedeb=ctreenew;
  nochange=1; movedist=0; 
  if(opt->SH_G>0 && opt->SH_G<nb-3) 
    maxmovedist=opt->SH_G;
  else
    maxmovedist=nb-3;

  if(print1)
    printf("\nStarting rearrangements\n");
  
  do{

    oldmovedist=movedist;

    if(nochange==1) movedist++;
    else movedist=1;

    if(print1 && !(nochange==0 && oldmovedist==1)){
      printf("\nCrossing %d internal branch", movedist);
      if(movedist>1) printf("es");
      printf("\n");
    }

    nochange=1;





    while(1){
      setclades(curtree, nb, seqname, l1, list1, l2, list2, dclade, gclade, ddist, gdist, nbdclade, nbgclade, dcladename, gcladename);
      restart_d=0;
      for(i=0;i<nbdclade;i++){
        for(j=0;j<nbdclade;j++){
	  if(abs(ddist[i][j])!=movedist) continue; 
	  if(ddist[i][j]<0) continue;
	  moveclade(dclade, newdclade, nbdclade, nb, ddist, i, j);
	  settree(newdclade, nbdclade, gclade, nbgclade, nb, newtree);

	  if(deja_evalue(newtree, nb, list_tree, nblist, nblistmax))
	    continue;
          if(solid_branch_crossed(newtree, curtree, nb, solid))
	    continue;

	  addtolist(newtree, nb, list_tree, nblist, nblistmax);
  	  nblist++;
	  if(print1){
	    printf("\nMoving %s toward %s\n", dcladename[i], dcladename[j]);
          }
	  if(print2)
	    printf("\n");
	  lkh=maxlike(nb, seq, seqname, newtree, NULL, NULL, nbbi, l1, list1, l2, list2, opt, ctreenew, NULL, NULL);

	  if(lkh>maxlkh){
  	    ctree_noblbs(ctreenew, ctreenew_nobl, strlen(ctreenew));
	    if(print1)
	      printf("New tree is optimal : \n%s\n\nRestarting rearrangements\n", ctreenew_nobl);
	    maxlkh=lkh;
	    copytree(newtree, curtree, nb);
            if(opt->SH_RESTART) save_current_best("current_best_tree", ctreenew_nobl, lkh, 0);
            if(opt->SH_RESTART) save_evaluated("evaluated_trees", ctreenew_nobl, lkh, 0);
	    restart_d=1;
	    nochange=0;
	    ctreenew=treedeb;
	    while(*ctreenew) {*ctreenew=0; ctreenew++; }
	    ctreenew=treedeb;
	    break;
	  }
	  else{
	    if(print1)
	      printf("No improvement:\n");
            ctree_noblbs(ctreenew, ctreenew_nobl, strlen(ctreenew));
            if(print1)
              printf("%s\n", ctreenew_nobl);
            if(opt->SH_RESTART) save_evaluated("evaluated_trees", ctreenew_nobl, lkh, 0);
	  }
	  ctreenew=treedeb;
	  while(*ctreenew) {*ctreenew=0; ctreenew++; }
	  ctreenew=treedeb;
        }
      if(restart_d) break;
      }
      if(!restart_d || movedist>1) break;
    }

    if(restart_d && movedist>1) continue;

    while(1){
      setclades(curtree, nb, seqname, l1, list1, l2, list2, dclade, gclade, ddist, gdist, nbdclade, nbgclade, dcladename, gcladename);
      restart_g=0;
      for(i=0;i<nbgclade;i++){
        for(j=0;j<nbgclade;j++){
	  if(abs(gdist[i][j])!=movedist) continue;
	  if(gdist[i][j]<0) continue;

	  moveclade(gclade, newgclade, nbgclade, nb, gdist, i, j);
	  settree(dclade, nbdclade, newgclade, nbgclade, nb, newtree);

	  if(deja_evalue(newtree, nb, list_tree, nblist, nblistmax))
	    continue;
          if(solid_branch_crossed(newtree, curtree, nb, solid))
	    continue;
	  addtolist(newtree, nb, list_tree, nblist, nblistmax);
	  nblist++;
	  if(print1){
	    printf("\nMoving %s toward %s\n", gcladename[i], gcladename[j]);
          }
	  if(print2){
	    printf("\n");
          }

	  lkh=maxlike(nb, seq, seqname, newtree, NULL, NULL, nbbi, l1, list1, l2, list2, opt, ctreenew, NULL, NULL);

	  if(lkh>maxlkh){
  	    ctree_noblbs(ctreenew, ctreenew_nobl, strlen(ctreenew));
	    if(print1)
	      printf("New tree is optimal : \n%s\n\nRestarting rearrangements\n", ctreenew_nobl);
	    maxlkh=lkh;
	    copytree(newtree, curtree, nb);
	    if(opt->SH_RESTART) save_current_best("current_best_tree", ctreenew_nobl, lkh, 0);
            if(opt->SH_RESTART) save_evaluated("evaluated_trees", ctreenew_nobl, lkh, 0);
            restart_g=1;
	    nochange=0;
	    ctreenew=treedeb;
	    while(*ctreenew) {*ctreenew=0; ctreenew++; }
	    ctreenew=treedeb;
	    break;
	  }
	  else{
	    if(print1)
	      printf("No improvement:\n");
            ctree_noblbs(ctreenew, ctreenew_nobl, strlen(ctreenew));
            if(print1)
              printf("%s\n", ctreenew_nobl);
            if(opt->SH_RESTART) save_evaluated("evaluated_trees", ctreenew_nobl, lkh, 0);
	  }
	  ctreenew=treedeb;
	  while(*ctreenew) {*ctreenew=0; ctreenew++; }
	  ctreenew=treedeb;
        }
      if(restart_g) break;
      }
      if(!restart_g || movedist>1) break;
    }
  } while(nochange==0 || movedist!=maxmovedist);



  	/* FINAL EVALUATION */

  if(print1)
    printf("\nFinal evaluation\n");
  if(print2)
    printf("\n");
  opt->print=trueprint;

  maxlkh=maxlike(nb, seq, seqname, curtree, NULL, NULL, nbbi, l1, list1, l2, list2, opt, NULL, ctree1, ctree2);

  if(ctree1){
    outfile1=fopen("treefile.eqgc", "w");
    outfile2=fopen("treefile.ndgc", "w");
    if(outfile1==NULL || outfile2==NULL){ printf("Cannot write tree file\n"); exit(0); }
    fprintf(outfile1, "%s\n", ctree1);
    fprintf(outfile2, "%s\n", ctree2);
    if(print1){
      printf("Tree is written into files : treefile.eqgc (equilibrium G+C content)\n");
      printf("                             treefile.ndgc (G+C content at each node)\n\n");
    }
  }


  free(lgbp); free(sortedlgbp); free(lgbi);
  free(bootvals);
  for(i=0;i<nb;i++){ free(nom[i]); free(ttree[i]); }
  free(list1); free(list2); free(curtree);

  
  return maxlkh;
}

 


main(int argc, char** argv){

/* 		A->0 ; C->1 ; G->2 ; T->3 		*/

  int i, j, nbseq, nb1, nb2, prov, secondes, minutes, heures, jours, muet, nbeval;
  char nomfinseq[50], nomfintree[50], nomfopt[50], *seq[MAXNSP], *seqname[MAXNSP], *comments[MAXNSP];
  char *ctree, *current_best, **evaluated, *ret;
  double maxl, runtime;
  FILE *treefile,*in, *optfile, *current_best_file, *evaluated_file;
  options opt;
  time_t debut, fin;

/* debugf=fopen("shake_debug", "w");*/


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
      exit(0);
    }
    sprintf(nomfinseq, "%s", argv[1]);
    muet=1;
  }
  fclose(in);

  if(muet)
    nbseq=readmasemuet(nomfinseq, seq, seqname, comments, MAXDATASET*MAXNSP);
  else
    nbseq=readmaseseqs(nomfinseq, seq, seqname, comments, MAXDATASET*MAXNSP);

  refresh(seq, nbseq, 0);

/* debug_nom=seqname; */

   

   /** input tree **/

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
      exit(0);
    }
  }
  
  ctree=readtree(treefile, nbseq);


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
      exit(0);
    }
  }

  getoptions(&opt, optfile);


  /** input saved work **/

  if(opt->SH_RESTART && argc==5){
    evaluated_file=fopen(argv[4], "r");

    if(!evaluated_file){
      printf("Cannot find evaluated file %s\n", argv[5]);
      exit(0);
    }

    evaluated=check_alloc(MAX_NBLISTMAX, sizeof(char*));
    nbeval=0; j=3*MAXLNAME*nbseq;
    while(1){
      evaluated[nbeval]=check_alloc(j, sizeof(char));
      ret=fgets(evaluated[nbeval], j, evaluated_file);
      if(ret==NULL) break;
      nbeval++;
    }
  }
  else{
    evaluated=NULL; nbeval=0;
  }

  time(&debut);


    /** shake **/

  maxl=shake(nbseq, seq, seqname, ctree, opt, evaluated, nbeval);


    /** output **/

  if(!opt->print->PRINT1)
    printf("%f\n", maxl);

  time(&fin);

  runtime=difftime(fin, debut);

  if(opt->print->PRINT1){
    jours=(int)(runtime/86400.);
    heures=(int)((runtime-86400*jours)/3600.);
    minutes=(int)((runtime-86400*jours-3600*heures)/60.);
    secondes=(int)(runtime-86400*jours-3600*heures-60*minutes);
    printf("(running time : ");
    if(jours) printf("%d days, ", jours);
    if(jours || heures) printf("%d hours, ", heures);
    if(jours || heures || minutes) printf("%d minutes, ", minutes);
    printf("%d seconds)\n\n", secondes);
  }
  
}



