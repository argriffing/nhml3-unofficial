#include "nhmlg.h"


tree* s_tree;
compute_option compute_o;

noeud root;
extern int gamma_nbcl;

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
      printf("Bad tree\n");
      exit(EXIT_FAILURE);
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
    exit(EXIT_FAILURE);
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

/*
for(i=0;i<nbotu;i++){
  for(j=0;j<nbbi;j++) printf("%d ", sorted[i][j]);
  printf("\n");
}

for(i=0;i<2*nbotu-2;i++)
  tabnd[i]->indice=i;
for(i=0;i<2*nbotu-2;i++){
  nd1=tabnd[i];
  printf("nd %d : %d %d %d\n", nd1->indice, nd1->v1?nd1->v1->indice:mo1, nd1->v2?nd1->v2->indice:mo1, nd1->v3?nd1->v3->indice:mo1);
}
*/



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
      exit(EXIT_FAILURE);
    }
  }
  s_tree->nbseq=nbseq;

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
    if(lgbp[0]>-0.5) opt->init->INIT_ROOT=root->l1/(root->l1+root->l2);
    else opt->init->INIT_ROOT=0.5;
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
    stoc(s_tree->node, nbseq, 1, ctree, -1);

  if(ctree1)
    stoc(s_tree->node, 1, nbseq, ctree1, 0);

  if(ctree2){
    set_node_deduced_comp(root);
    stoc(s_tree->node, 1, nbseq, ctree2, 1);
  }


  /* FREEING */

  free(drapeau); free(drapeau2); free(site);
  free_tree(s_tree, gamma_nbcl);
  free(s_tree);
  free_node(root, lgseq, nbparam, 0, gamma_nbcl);

  return maxlkh;

}




int samebranch(int** sorted, int nb, int br1, int br2){

  int i, eq=0, diff=0;

  for(i=0;i<nb;i++){
    if(sorted[i][br1]==sorted[i][br2]) eq=1;
    else diff=1;
    if(diff&&eq) return 0;
  }
  return 1;
}



int compatible(int** sorted, int nb, int br1, int br2){

  int i, zu=0, zz=0, uu=0, uz=0;

  for(i=0;i<nb;i++){
    if(sorted[i][br1]==0 && sorted[i][br2]==0) zz=1;
    else if(sorted[i][br1]==0 && sorted[i][br2]==1) zu=1;
    else if(sorted[i][br1]==1 && sorted[i][br2]==0) uz=1;
    else if(sorted[i][br1]==1 && sorted[i][br2]==1) uu=1;
    if(zz && zu && uz && uu) return 0;
  }
  return 1;
}


int goodnewbranch(int** sorted, int nb, int nbbi){
  int l;

  for(l=0;l<nbbi-1;l++){
    if(!compatible(sorted, nb, l, nbbi-1))
      break;
  }
  if(l!=nbbi-1) return 0;

  for(l=0;l<nbbi-1;l++){
    if(samebranch(sorted, nb, l, nbbi-1))
      break;
  }
  if(l!=nbbi-1) return 0;

  return 1;

}


void setnewbranch(int** tree, int nbotu, int numbi, int* listotu, int* listotu_init, int listinitsize, int j, int k){
  int l;

  for(l=0;l<nbotu;l++) tree[l][numbi]=0;
  for(l=0;l<listinitsize;l++)
    if(listotu[l]==listotu[j] || listotu[l]==-listotu[j]-1 
      || listotu[l]==listotu[k] || listotu[l]==-listotu[k]-1)
      tree[listotu_init[l]][numbi]=1;

}



double split(int nb, char** seq, char** seqname, char* ctree, options opt){

  int nb2, i, j, k, l, *ttree[MAXNSP], **sorted, nbbi, l1, l2, nbcycle, nbotud, nbotug, *otud, *otug, otu1, otu2, *otud_init, *otug_init, maxnbpair, pres_grossiere=0, truepres, nbpair, nbbest, best, nbgrouping,
print1, print2;
  double *lgbp, *sortedlgbp, *lgbi, *bootvals, fracroot1, *lkh, maxlkh, inf_bound, prov;
  char *nom[MAXNSP], *groupnom[MAXNSP], **list1, **list2, racine, *whichpart, whichpart_final, *ctreenew, *ctree1, *ctree2, *provc;
  FILE* outfile1, *outfile2;
  typedef int pair[2];
  pair *tried;
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


		/* READ TREE STRING */

  lgbp=(double*)check_alloc(nb+1, sizeof(double));
  sortedlgbp=(double*)check_alloc(nb+1, sizeof(double));
  lgbi=(double*)check_alloc(nb+1, sizeof(double));
  for(i=0;i<nb+1;i++)
    lgbp[i]=lgbi[i]=-1.;
  bootvals=(double*)check_alloc(nb+1, sizeof(double));
  for(i=0;i<nb+1;i++){
    nom[i]=(char*)check_alloc(MAXLNAME+1, sizeof(char));
    groupnom[i]=(char*)check_alloc(nb*(MAXLNAME+3)+1, sizeof(char));
    ttree[i]=(int*)check_alloc(nb, sizeof(int));
  }
  list1=check_alloc(nb, sizeof(char*));
  list2=check_alloc(nb, sizeof(char*));
  sorted=(int**)check_alloc(nb+1, sizeof(int*));
  ctreenew=(char*)check_alloc(2*MAXLNAME*nb, sizeof(char));
  ctree1=(char*)check_alloc(2*MAXLNAME*nb, sizeof(char));
  ctree2=(char*)check_alloc(2*MAXLNAME*nb, sizeof(char));
  maxnbpair=nb*(nb+1)/2;
  tried=(pair*)check_alloc(maxnbpair, sizeof(pair));
  whichpart=(char*)check_alloc(maxnbpair, sizeof(char));
  lkh=(double*)check_alloc(maxnbpair, sizeof(double));

  for(i=0;i<nb;i++) sprintf(groupnom[i], "%s", seqname[i]);

  nb2=ctot(ctree, ttree, lgbi, lgbp, bootvals, nom, &racine, &nbbi);

  if(nb2<nb){
    printf("More species in sequence file than in tree file\n");
    exit(EXIT_FAILURE);
  }
  if(nb2>nb){
    printf("More species in tree file than in sequence file\n");
    exit(EXIT_FAILURE);
  }

  if(racine=='r'){
    unroot(ttree, nb, lgbi, lgbp, bootvals, nom, list1, list2, &l1, &l2, &fracroot1);
  }
  else{
    printf("Tree must be rooted\n");
    exit(EXIT_FAILURE);
  }
  nbbi--;

  for(i=0;i<nb;i++){
    for(j=0;j<nb;j++){
      if(samename(seqname[i], nom[j])){
	sorted[i]=ttree[j];
	sortedlgbp[i]=lgbp[j];
	break;
      }
    }
  }

  truepres=opt->converge->PRECISION;
  opt->converge->PRECISION=pres_grossiere;
  if(print1)
    printf("\nInitial evaluation\n");
  if(print2)
    printf("\n");

  maxlkh=maxlike(nb, seq, seqname, sorted, lgbi, sortedlgbp, nbbi, l1, list1, l2, list2, opt, NULL, NULL, NULL);

  if(print1)
    printf("\n");


		/* SET OTU LISTS */

  otud=(int*)check_alloc(l1, sizeof(int));
  otud_init=(int*)check_alloc(l1, sizeof(int));
  otug=(int*)check_alloc(l2, sizeof(int));
  otug_init=(int*)check_alloc(l2, sizeof(int));
  k=0;
  for(i=0;i<l1;i++){
    for(j=0;j<nb;j++){
      if(samename(list1[i], seqname[j])){
	otud[k]=j;
	k++;
	break;
      }
    }
  }
  if(k!=l1) printf("erreur otud\n");
  k=0;
  for(i=0;i<l2;i++){
    for(j=0;j<nb;j++){
      if(samename(list2[i], seqname[j])){
	otug[k]=j;
	k++;
	break;
      }
    }
  }
  if(k!=l2) printf("erreur otug\n");
  nbotud=l1; nbotug=l2;
  for(i=0;i<l1;i++) otud_init[i]=otud[i];
  for(i=0;i<l2;i++) otug_init[i]=otug[i];



		/* STAR DECOMPOSITION */

  nbcycle=nb-3-nbbi;
  for(i=0;i<nbcycle;i++){
    nbgrouping=nbotud*(nbotud-1)/2 + nbotug*(nbotug-1)/2;
    if(nbotud==2) nbgrouping--;
    if(nbotug==2) nbgrouping--;
    if(print1)
      printf("STEP %d : (%d groupings to try)\n\n", i+1, nbgrouping);
    nbbi++; otu1=otu2=-1;
    nbpair=0;
    opt->converge->PRECISION=pres_grossiere;

    /* groupements d'otu droits (evaluation grossiere) */
    if(nbotud>2){
      for(j=0;j<l1;j++){
        if(otud[j]<0) continue;
        for(k=j+1;k<l1;k++){
	  if(otud[k]<0) continue;
	  setnewbranch(sorted, nb, nbbi-1, otud, otud_init, l1, j, k);
	  if(!goodnewbranch(sorted, nb, nbbi)) continue; 
          if(print1) 
    	    printf("grouping %s and %s \n", groupnom[otud[j]], groupnom[otud[k]]);
	  if(print2)
	    printf("\n");
	  lkh[nbpair]=maxlike(nb, seq, seqname, sorted, NULL, NULL, nbbi, l1, list1, l2, list2, opt, NULL, NULL, NULL);
	  tried[nbpair][0]=j;
	  tried[nbpair][1]=k;
	  whichpart[nbpair]='d';
	  nbpair++;
        }
      }
    }

    /* groupements d'otu gauches */
    if(nbotug>2){
      for(j=0;j<l2;j++){
        if(otug[j]<0) continue;
        for(k=j+1;k<l2;k++){
	  if(otug[k]<0) continue;
	  setnewbranch(sorted, nb, nbbi-1, otug, otug_init, l2, j, k);
	  if(!goodnewbranch(sorted, nb, nbbi)) continue; 
          if(print1)
	    printf("grouping %s and %s \n", groupnom[otug[j]], groupnom[otug[k]]);
	  if(print2)
	    printf("\n");
	  lkh[nbpair]=maxlike(nb, seq, seqname, sorted, NULL, NULL, nbbi, l1, list1, l2, list2, opt, NULL, NULL, NULL);
	  tried[nbpair][0]=j;
	  tried[nbpair][1]=k;
	  whichpart[nbpair]='g';
	  nbpair++;
	}
      }
    }



    /* departage fin des meilleurs */

    if(nbpair==0){
      printf("Error in split decomposition\n");
      exit(EXIT_FAILURE);
    }
    opt->converge->PRECISION=truepres;
    maxlkh=lkh[0];
    for(l=1;l<nbpair;l++)
      if(lkh[l]>maxlkh) maxlkh=lkh[l];
    inf_bound=maxlkh-pow(10., -(double)pres_grossiere);
    nbbest=0;
    for(l=0;l<nbpair;l++){
      if(lkh[l]>=inf_bound){
	nbbest++;
        best=l;
     	whichpart_final=whichpart[l];
      }
    }

    if(whichpart_final=='d'){otu1=otud[tried[best][0]]; otu2=otud[tried[best][1]];}
    else {otu1=otug[tried[best][0]]; otu2=otug[tried[best][1]];}

    if(nbbest>1){
      if(print1)
	printf("Precise evaluation of %d tied alternative groupings\n\n", nbbest);
       k=0;
      for(l=0;l<nbpair;l++){
        if(lkh[l]>=inf_bound){
	  if(print1) printf("tied %d:\n", k+1);
          if(print2)
	    printf("\n");
	  k++;
	  if(whichpart[l]=='d')
	    setnewbranch(sorted, nb, nbbi-1, otud, otud_init, l1, tried[l][0], tried[l][1]);
	  else
	    setnewbranch(sorted, nb, nbbi-1, otug, otug_init, l2, tried[l][0], tried[l][1]);
	  prov=maxlike(nb, seq, seqname, sorted, NULL, NULL, nbbi, l1, list1, l2, list2, opt, NULL, NULL, NULL);
	  if(prov>maxlkh){
	    maxlkh=prov;
	    whichpart_final=whichpart[l];
	    if(whichpart[l]=='d'){otu1=otud[tried[l][0]]; otu2=otud[tried[l][1]];}
	    else {otu1=otug[tried[l][0]]; otu2=otug[tried[l][1]];}
 	  }
        }
      }
    }
    
    

    /* modification de l'arbre, mise a jour des listes */

    if(nbbest==0){
      printf("No pairwise grouping increases the likelihood: program error\n");
      return maxlkh;
    }
    else {
      if(print1) 
        printf("Best grouping for step %d:\n %s and %s\n\n", i+1, groupnom[otu1], groupnom[otu2]);
    }


    if(whichpart_final=='d'){
      nbotud--;
      for(l=0;l<nb;l++)
	sorted[l][nbbi-1]=0;
      for(l=0;l<l1;l++){
	if(otud[l]==otu1 || otud[l]==otu2 || otud[l]==-otu1-1 || otud[l]==-otu2-1)
	  sorted[otud_init[l]][nbbi-1]=1;
      }
      provc=groupnom[otu1];
      for(l=0;l<l1;l++){
	if(otud[l]==otu2 || otud[l]==-otu2-1){
	  otud[l]=-otu1-1;
	  while(*provc) provc++;
	  sprintf(provc, "+%s", seqname[otud_init[l]]);
	}
      }
    }
    if(whichpart_final=='g'){
      nbotug--;
      for(l=0;l<nb;l++)
	sorted[l][nbbi-1]=0;
      provc=groupnom[otu1];
      for(l=0;l<l2;l++){
	if(otug[l]==otu1 || otug[l]==otu2 || otug[l]==-otu1-1 || otug[l]==-otu2-1)
	  sorted[otug_init[l]][nbbi-1]=1;
      }
      for(l=0;l<l2;l++){
	if(otug[l]==otu2 || otug[l]==-otu2-1){
	  otug[l]=-otu1-1;
	  while(*provc) provc++;
	  sprintf(provc, "+%s", seqname[otug_init[l]]);
	}
      }
    }
  }

  	/* FINAL EVALUATION */

  if(print1)
    printf("\nFinal evaluation\n");
  if(print2)
    printf("\n");
  opt->converge->PRECISION=truepres;
  opt->print=trueprint;
  maxlkh=maxlike(nb, seq, seqname, sorted, NULL, NULL, nbbi, l1, list1, l2, list2, opt, NULL, ctree1, ctree2);



  if(ctree1){
    outfile1=fopen("treefile.eqgc", "w");
    outfile2=fopen("treefile.ndgc", "w");
    if(outfile1==NULL || outfile2==NULL){ printf("Cannot write tree file\n"); exit(EXIT_FAILURE); }
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
  free(list1); free(list2); free(sorted);
  free(otud); free(otug); free(otud_init); free(otug_init);


  return maxlkh;
}




int main(int argc, char** argv){

/* 		A->0 ; C->1 ; G->2 ; T->3 		*/

  int i, j, nbseq, nb1, nb2, prov, secondes, minutes, heures, jours, muet;
  char nomfinseq[50], nomfintree[50], nomfopt[50], *seq[MAXNSP], *seqname[MAXNSP], *comments[MAXNSP];
  char *ctree;
  double maxl, runtime;
  FILE *treefile,*in, *optfile;
  options opt;
  time_t debut, fin;

  time(&debut);

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

  refresh(seq, nbseq, 0);
   

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
      exit(EXIT_FAILURE);
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
      exit(EXIT_FAILURE);
    }
  }

  getoptions(&opt, optfile);



    /** star decomposition **/

  maxl=split(nbseq, seq, seqname, ctree, opt);


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

  return 0;
}
