/***************************************************************************/
/*                                                                         */  
/*                                                                         */  
/*                               ANCHOR                                    */
/*                                                                         */  
/*      Prediction of Protein Binding Regions in Disordered Proteins       */  
/*         Balint Meszaros, Istvan Simon and Zsuzsanna Dosztanyi           */
/*          Institute of Enzymology, Biological Research Center            */
/*           Hungarian Academy of Sciences, Budapest, Hungary              */  
/*                                                                         */
/***************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#define max_line_len 1000
#define AA "GAVLIFPSTCMWYNQDEKRH"
#define AAN 20
#define ML 1000
#define MSL 10000
#define WINDOW2 7
#define WINDOW_IUPRED 100
#define max_reg 1000
#define pi 3.141592654

/* Optimized parameters:*/

int W_COMP_1=60;
int W_IUPR_1=25;
double E_COEFF_1=-3.847;
double ED_COEFF_1=7.9853;
double IUP_COEFF_1=4.6297;

int W_COMP_2=60;
int W_IUPR_2=27;
double E_COEFF_2=-4.149;
double ED_COEFF_2=6.7728;
double IUP_COEFF_2=6.0754;

int W_COMP_3=90;
int W_IUPR_3=29;
double E_COEFF_3=-4.585;
double ED_COEFF_3=5.4879;
double IUP_COEFF_3=6.9896;


double CUTOFF_1=5.930970;
double CUTOFF_2=6.630499;
double CUTOFF_3=7.093414;


typedef struct
{
	char name[1000];
	char *seq;
	int len;
} SEQ_STR;

typedef struct
{
	double *global;
	double *local;
	double *external;
	double *locglob;
	double *extglob;
	double *glob_smooth;
	double *iupred;
	double *crit_e;
	double *crit_ed;
	int *crit_iup;
	double *crit_iup_bal;
	double *crit_iup_jobb;
	double *iup_av;
	int *prediction;
	double *score;
	int *score2;
} E_STR;

typedef struct 
{
  double *distro;
  double min, max;
  double step;
  double cutoff;
  int nb;
  double *comp_ext;  
  double **E_mat;
} P_STR;

P_STR *Read_Data(char *path, char *fn);
void Usage(char *prgname);
E_STR *IUPred(SEQ_STR *SEQ,P_STR *P,int WI);
SEQ_STR  *Get_Seq(char *name);
void Smooth(double *data,int length,int win_size);
void Combine(E_STR *E,int l,double E_COEFF,double ED_COEFF,double IUP_COEFF,double CUTOFF,int W);
int Get_Regions(double *data,int length,int **r);
void Check_Regions(int **r,int numreg,double *iupred);
double convert_to_p_value(double score);

void *Gr_malloc(size_t size);
double **DMatrix(int n_rows, int n_cols);

int main(int argc, char **argv)
{
  P_STR *P;
  SEQ_STR *SEQ;
  E_STR *E1;
  E_STR *E2;
  E_STR *E3;
  E_STR *E_IUPred;
  int i,j,num_reg,**regions;
  double *p;
  int *reg_filt;
  char *path;
  int verb=0;
  char sp_name[1000];
  int nrt,nrf;

  if (argc<2) {Usage(argv[0]);}
  path=NULL;  

  if (argc>=2) {
    for (i=1;i<argc;i++) {
      if (strncmp(argv[i],"-v",2)==0) verb=1;
      else if (strncmp(argv[i],"-d",2)==0) {
        i++;
        if (i>=argc) printf("Data directory is missing\n"),exit(1);
        path=calloc(strlen(argv[i]),sizeof(char));
        sprintf(path,"%s",argv[i]);
      }
      else sprintf(sp_name,"%s",argv[i]);             
    }
  }
  if (path==NULL) {
    if ((path=getenv("ANCHOR_PATH"))==NULL) {
      path=calloc(200,sizeof(char));
      sprintf(path,"%s","./");
    } 
  }


  if (path==NULL) printf("No datadir was specified\n"),exit(1);


  

  P=Read_Data(path,"anchordata");

        
  SEQ=Get_Seq(sp_name);





  E_IUPred=IUPred(SEQ,P,WINDOW_IUPRED);
  E1=IUPred(SEQ,P,W_COMP_1);
  E2=IUPred(SEQ,P,W_COMP_2);
  E3=IUPred(SEQ,P,W_COMP_3);




  Smooth(E1->local,SEQ->len,4);
  Smooth(E1->local,SEQ->len,4);
  Smooth(E1->external,SEQ->len,4);
  Smooth(E1->external,SEQ->len,4);
  Smooth(E1->locglob,SEQ->len,4);
  Smooth(E1->locglob,SEQ->len,4);
  Smooth(E1->extglob,SEQ->len,4);
  Smooth(E1->extglob,SEQ->len,4);
  
  Smooth(E2->local,SEQ->len,4);
  Smooth(E2->local,SEQ->len,4);
  Smooth(E2->external,SEQ->len,4);
  Smooth(E2->external,SEQ->len,4);
  Smooth(E2->locglob,SEQ->len,4);
  Smooth(E2->locglob,SEQ->len,4);
  Smooth(E2->extglob,SEQ->len,4);
  Smooth(E2->extglob,SEQ->len,4);
  
  Smooth(E3->local,SEQ->len,4);
  Smooth(E3->local,SEQ->len,4);
  Smooth(E3->external,SEQ->len,4);
  Smooth(E3->external,SEQ->len,4);
  Smooth(E3->locglob,SEQ->len,4);
  Smooth(E3->locglob,SEQ->len,4);
  Smooth(E3->extglob,SEQ->len,4);
  Smooth(E3->extglob,SEQ->len,4);
  

  Combine(E1,SEQ->len,E_COEFF_1,ED_COEFF_1,IUP_COEFF_1,CUTOFF_1,W_IUPR_1);
  Combine(E2,SEQ->len,E_COEFF_2,ED_COEFF_2,IUP_COEFF_2,CUTOFF_2,W_IUPR_2);
  Combine(E3,SEQ->len,E_COEFF_3,ED_COEFF_3,IUP_COEFF_3,CUTOFF_3,W_IUPR_3);
  

  p=calloc(SEQ->len,sizeof(double));
  reg_filt=calloc(SEQ->len,sizeof(int));

  for (i=0;i<SEQ->len;i++) {
    p[i]=convert_to_p_value((E1->score[i]+E2->score[i]+E3->score[i])/3);
  }

  regions=malloc(4*sizeof(int *));
  for(i=0;i<4;i++) regions[i]=malloc(max_reg*sizeof(int));
  num_reg=Get_Regions(p,SEQ->len,regions);
  Check_Regions(regions,num_reg,E_IUPred->iupred);


  for (i=0;i<SEQ->len;i++) {reg_filt[i]=0;}

  for (i=0;i<SEQ->len;i++) {
    for(j=0;j<num_reg;j++) {
      if((i>=regions[0][j]-1)&&(i<regions[1][j])&&(regions[3][j]==0)) 
	reg_filt[i]=1;
    }
  }


  
  
  printf("# ANCHOR \n");
  printf("# Copyright (c)  ");
  printf("B. Meszaros, I. Simon and Z. Dosztanyi, 2009\n");
  printf("# PLoS Comput Biol 5(5): e1000376. doi:10.1371/journal.pcbi.1000376\n");
  printf("#\n");

  nrt=0;nrf=0;
  for(i=0;i<num_reg;i++) {
    if (regions[3][i]==0)   nrt++;
    else nrf++;    
  }

  printf("# Predicted binding regions\n");
  if(nrt==0) {printf("#   None.\n");}
  else {
    printf("#\t%-15s\t%8s\t%8s\t%8s\n","No.","Start","End","Length");
  }
 
  for(i=0,j=0;i<num_reg;i++) {
    if (regions[3][i]==0) {
      printf("#\t%-15d\t%8d\t%8d\t%8d\n",j+1,regions[0][i],
	     regions[1][i],regions[2][i]);
      j++;
    }
  }

  if (nrf>0) { 
    printf("#\n# Filtered regions\n");
    printf("#\t%-15s\t%8s\t%8s\t%8s\n","No.","Start","End","Length");
    
    for(i=0,j=0;i<num_reg;i++) {
      if (regions[3][i]!=0) {
	printf("#\t%-15d\t%8d\t%8d\t%8d\n",j+1,regions[0][i],
	       regions[1][i],regions[2][i]);
	j++;
      }
    }
  }


  if (verb==0) {
    printf("#\n# Prediction profile output:\n");
    printf("#   Columns:\n#   1 - Amino acid number\n#   2 - One letter code\n");
    printf("#   3 - ANCHOR probability value\n");
    printf("#   4 - ANCHOR output\n");
    printf("#\n");
    for (i=0;i<SEQ->len;i++) {
      printf("%d\t%c\t%8.4f\t%5d\n",i+1,SEQ->seq[i],p[i],reg_filt[i]);
    }
  }

  if (verb==1) {
    printf("#\n# Prediction profile output:\n");
    printf("#   Columns:\n#   1 - Amino acid number\n");
    printf("#   2 - One letter code\n");
    printf("#   3 - ANCHOR probability value\n");
    printf("#   4 - ANCHOR output\n");
    printf("#   5 - IUPred probability value\n");
    printf("#   6 - ANCHOR score \n");
    printf("#   7 - S \n");
    printf("#   8 - Eint \n");
    printf("#   9 - Egain \n");

    for (i=0;i<SEQ->len;i++) {
      printf("%5d\t%c\t%9.4f\t%5d\t%9.4f\t%9.4f\t%9.4f\t%9.4f\t%9.4f\n",
	     i+1,SEQ->seq[i],
             p[i],
	     reg_filt[i],
	     E_IUPred->iupred[i],
	     (E1->score[i]+E2->score[i]+E3->score[i])/3,
	     (E1->iup_av[i]+E2->iup_av[i]+E3->iup_av[i])/3,
	     -(E1->external[i]+E2->external[i]+E3->external[i])/3,
           (E1->extglob[i]+E2->extglob[i]+E3->extglob[i])/3);
    }
  }




  exit(0);
}


double convert_to_p_value(double score)
{
  double mu_pos=0.830866,mu_neg=-2.835816;
  double sigma_pos=2.387250,sigma_neg=1.946386;
  double neg_factor=2.21807;
  double p,val_pos,val_neg;
  
  if (score<-90) return 0;

  val_pos=(1.0/(sigma_pos*(sqrt(2*pi))))*exp(-(score-mu_pos)*(score-mu_pos)/(2*sigma_pos*sigma_pos));
  val_neg=(1.0/(sigma_neg*(sqrt(2*pi))))*exp(-(score-mu_neg)*(score-mu_neg)/(2*sigma_neg*sigma_neg));

  
  p=1.0/(1.0+((neg_factor*val_neg)/val_pos));
  if(score<(mu_pos-3*sigma_pos)) {p*=exp(score-(mu_pos-3*sigma_pos));}
  
  return(p);
}


void Combine(E_STR *E,int l,double E_COEFF,double ED_COEFF,double IUP_COEFF,double CUTOFF,int W)
{
  int i,j;
  double jobb,bal,j_b;
  int jobb_num,bal_num;
  
  E->iup_av=calloc(l,sizeof(double));
  E->prediction=calloc(l,sizeof(int));
  E->score=calloc(l,sizeof(double));
  
  for (i=0;i<l;i++) {
    jobb=0,bal=0;
    jobb_num=0,bal_num=0;
    for (j=i-W;j<i+W;j++) {
      if ((j>=0)&&(j<i)) {
	bal+=E->iupred[j];
	bal_num++;
      } 
      if ((j>i)&&(j<l)) {
	jobb+=E->iupred[j];
	jobb_num++;
      }
    }
    j_b=(bal+jobb)/(bal_num+jobb_num);
    
    E->iup_av[i]=j_b;
    
    E->score[i]=E->external[i]*E_COEFF+E->extglob[i]*ED_COEFF+E->iup_av[i]*IUP_COEFF-CUTOFF;
    /*		printf("%d\t%f\n",i,E->score[i]);*/
  }
}


int Get_Regions(double *data,int length,int **r)
{
  int i;
  int num_reg=0,on=0;
  
  for (i=0;i<length;i++) {
    if ((i==0)&&(data[i]>0.5)) {
      num_reg=1;
      on=1;
      r[0][num_reg-1]=1;
    }
    
    if ((i==length-1)&&(on==1)) {
      r[1][num_reg-1]=length;
      r[2][num_reg-1]=r[1][num_reg-1]-r[0][num_reg-1]+1;
    }
    
    if ((0<i)&&(i<length-1)) {
      if ((data[i]>0.5)&&(data[i-1]<=0.5)) {
	num_reg++;
	r[0][num_reg-1]=i+1;
	on=1;
			}
      if ((data[i]<=0.5)&&(data[i-1]>0.5)) {
	r[1][num_reg-1]=i;
	r[2][num_reg-1]=r[1][num_reg-1]-r[0][num_reg-1]+1;
	on=0;
      }
    }
  }
  return(num_reg);	

}

void Check_Regions(int **r,int numreg,double *iupred)
{
  int i,j;
  double min_score;

  for(i=0;i<numreg;i++) {
    min_score=1;
    r[3][i]=0;

    for(j=r[0][i]-1;j<r[1][i];j++) {
      if(iupred[j]<min_score) min_score=iupred[j];
    }


    if(r[2][i]<6) r[3][i]=1;
    if((min_score<0.09)&&(r[2][i]>=6)) r[3][i]=2;
  }
}


void Smooth(double *data,int length,int win_size)
{
  int i,j,begin,end;
  double value,*output;
  int points;
  
  output=calloc(length,sizeof(double));
  
  for (i=0;i<length;i++) {
    if (i>win_size) begin=i-win_size;
    else begin=0;
    if (i<length-1-win_size) end=i+win_size;
    else end=length-1;
    points=0;
    value=0;
    for (j=begin;j<=end;j++) {value+=data[j];points++;}
    value/=points;
    if (i>=length) printf("%5d\n",i);
    output[i]=value;
  }
  for (i=0;i<length;i++) data[i]=output[i];

  free(output);
}


E_STR *IUPred(SEQ_STR *SEQ,P_STR *P,int WINDOW )
{
  int i,j,k;
  int ind;
  int p;
  int upper,lower;
  int *comp,*comp_local;
  int n_glob=0,n_loc=0;
  double min,max,step;
  E_STR *E;

  min=P->min;
  max=P->max;
  step=P->step;

  E=malloc(sizeof(E_STR));  
  E->global=calloc(SEQ->len,sizeof(double));
  E->glob_smooth=calloc(SEQ->len,sizeof(double));
  E->iupred=calloc(SEQ->len,sizeof(double));
  E->local=calloc(SEQ->len,sizeof(double));
  E->external=calloc(SEQ->len,sizeof(double));
  E->locglob=calloc(SEQ->len,sizeof(double));
  E->extglob=calloc(SEQ->len,sizeof(double));
  
  comp=calloc(AAN,sizeof(int));
  comp_local=calloc(AAN,sizeof(int));

   
  for (j=0;j<SEQ->len;j++) {
    
    n_glob=0;
    n_loc=0;
    if (j<WINDOW+1) {lower=0;}
    else {lower=j-WINDOW;}
    if ((SEQ->len-j)<WINDOW+1) {upper=SEQ->len-1;}
    else {upper=j+WINDOW;}
    for (k=0;k<AAN;k++) {comp[k]=0,comp_local[k]=0;}
    for (k=lower;k<upper+1;k++) {
      if ((k!=j)&&(k!=j-1)&&(k!=j+1)) {
	ind=strchr(AA,toupper(SEQ->seq[k]))-AA;
	if ((ind>-1)&&(ind<AAN)) {
	  comp[ind]++;
	  n_glob++;
	  if (((k-j)>=(-WINDOW2))&&((k-j)<=(WINDOW2))) {
	    comp_local[ind]++;
	    n_loc++;
	  }
	}
      }
    }
    
    ind=strchr(AA,toupper(SEQ->seq[j]))-AA;
    if ((ind<0)||(ind>=AAN)) continue;

    for (k=0;k<AAN;k++) {
      E->global[j]+=(P->E_mat[ind][k]*comp[k]);
      E->local[j]+=(P->E_mat[ind][k]*comp_local[k]);
      E->external[j]+=(P->E_mat[ind][k]*P->comp_ext[k]);
    }
		
    /*printf("%d %d %f\n",j+1,n_glob,E->global[j]);
     */
    E->global[j]/=n_glob;
    E->glob_smooth[j]=E->global[j];
    E->local[j]/=n_loc;
    
    E->locglob[j]=E->global[j]-E->local[j];
    E->extglob[j]=E->global[j]-E->external[j];
  }
  
  Smooth(E->glob_smooth,SEQ->len,10);

  for (i=0;i<SEQ->len;i++) {
    if (-E->glob_smooth[i]<=min+2*step) E->iupred[i]=1;
    if (-E->glob_smooth[i]>=max-2*step) E->iupred[i]=0;
    if ((-E->glob_smooth[i]>min+2*step)&&(-E->glob_smooth[i]<max-2*step)) {
      p=(int)((-E->glob_smooth[i]-min)*(1.0/step));     
      E->iupred[i]=P->distro[p];
    }
  }

  return E;

}




SEQ_STR *Get_Seq(char *fn)
{
  char line[ML];
  char c=0;
  int j;
  FILE *f;
  SEQ_STR *SEQ;



  if ((fn==NULL)||(strlen(fn)==0)) {
    printf("No sequence filename\n"),exit(1);
  }

  if ((f=fopen(fn,"r"))==NULL) {
    printf("Could not open %s\n",fn),exit(1);
  }

  SEQ=malloc(sizeof(SEQ_STR));
  SEQ->seq=calloc(MSL,sizeof(char));

  if (fgets(line,ML,f)==NULL) exit(1);
  /*sscanf(&line[1],"%s",SEQ->name);*/
  sprintf(SEQ->name,"%s",fn);

  j=0;
  while ((c!='>') && (!feof(f))) {
    c=fgetc(f);
    if (isalpha((int)(c))) {
      SEQ->seq[j]=c;
      j++;
      if ((j>0)&&(j%MSL==0)) {
        SEQ->seq=realloc(SEQ->seq,(j+MSL)*sizeof(char));
      }      
    }    
  }
  SEQ->seq[j]=0;
  SEQ->len=j;

  
#ifdef DEBUG
  printf("%s %5d\n%s\n",SEQ->name,SEQ->len,SEQ->seq);
#endif
  


  return SEQ;
}




P_STR *Read_Data(char *path, char *fn)
{
  FILE *f;
  char ln[max_line_len];
  int i,nb,set;
  double v,min,max,cutoff,c; 
  P_STR *P;
  char *fullfn;
  int sl;
  int id1,id2;
  double val;
  double **E_mat;
  double value;
  char s[10];
  int p,ind;
  double *com;

  sl=strlen(path)+strlen(fn)+2;
  fullfn=malloc(sl*sizeof(char));
  sprintf(fullfn,"%s/%s",path,fn);


  if ((f=fopen(fullfn,"r"))==NULL) {
    printf("Could not open %s\n",fullfn);
    exit(1);
  }

  P=malloc(sizeof(P_STR));


  /* Read energy predictor matrix  */
  E_mat=DMatrix(AAN,AAN);
  for (i=0;i<AAN*AAN;i++) {
    if (fgets(ln,max_line_len,f)==NULL) break;
    sscanf(ln,"%d%*s%d%*s%lf",&id1,&id2,&val);
    /*printf("%5d %5d %10.4f\n",id1,id2,val);*/
    E_mat[id1][id2]=val;
  }
  P->E_mat=E_mat;
  /* Read composition   */ 

  com=calloc(AAN,sizeof(double));
  for (i=0;i<AAN;i++) {
    if (fgets(ln,max_line_len,f)==NULL) break;
    sscanf(ln,"%d%s%lf\n",&p,s,&value);  
    /*printf("%d %s %f\n",p,s,value);*/
    ind=strchr(AA,s[0])-AA;
    com[ind]=value;
  }
  P->comp_ext=com;


  /* Read histogram */
  if (fscanf(f,"%*s %lf %lf %d\n",&min, &max, &nb)==0) exit(1);
  
  P->distro=malloc(nb*sizeof(double ));
  for (i=0;i<nb;i++) P->distro[i]=0;
  
  cutoff=-1000;
  for (i=0,set=0;i<nb;i++) {
    if (fgets(ln,max_line_len,f)==0) exit(1);
    if (feof(f)) break;
    if (ln[0]=='#') continue;
    if (sscanf(ln,"%*s %lf %*s %*s   %lf\n", &c,&v)==0) exit(1);
    if ((set==0)&&(v<=0.5)) {set=1;cutoff=c;}
    P->distro[i]=v;
  }

  fclose(f);
	
  if (cutoff==-1000) {printf("Cutoff is not set\n"),exit(1);}

  P->max=max;
  P->min=min;
  P->nb=nb;
  P->cutoff=cutoff;	
  
  P->step=(max-min)/nb;
  P->cutoff-=P->step;

  free(fullfn);
  
  return P;

}

void Usage(char *prgname)
{
  printf("Usage: %s sequence file [-v] [-d data_dir]\n",prgname);
  printf("  sequence file  - filename containing a single fasta formatted sequence\n");
  printf("  options\n");
  printf("  -v          : produce detailed profile output\n");
  printf("  -d data dir : specify the location of directory containing anchordata file\n");

  exit(-1);
}








double **DMatrix(int n_rows, int n_cols)
{
    double **matrix;
    int   i;
 
    matrix = (double **) Gr_malloc(n_rows*sizeof(double *));
    matrix[0] = (double *) Gr_malloc(n_rows*n_cols*sizeof(double));
 
    for (i = 1; i < n_rows; i++)
      matrix[i] = matrix[i-1] + n_cols;
    return matrix;

}


void *Gr_malloc(size_t size)
{
 
    void    *new_mem;
    
    if (size == 0)
      return NULL;
    new_mem = malloc(size);
    if (new_mem == NULL) {
      fprintf(stderr, "can't allocate enough memory: %zu bytes\n", size);
    }
    return new_mem;                                                             
}
