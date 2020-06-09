#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<math.h>
#include<time.h>

#define FILENAME_INPUT "in.bin"
#define FILENAME_OUTPUT "out_dqds_nakabayashi.bin"
#define MAT(A,I,J,LDA) (A[(LDA)*(J)+(I)])
#define MAX(a,b) ((a)>=(b) ? (a):(b))
#define MIN(a,b) ((a)<=(b) ? (a):(b))

int has_non_positive(int n, double *x)
{
  int i;
  for(i=0; i<n; i++){
    if(x[i]<=0){ return 1; }
  }
  return 0;
}


int has_nan_or_inf(int n, double *x)
{
  int i;
  for(i=0; i<n; i++){
    if(!isfinite(x[i])){ return 1; }
  }
  return 0;
}

void print_vector(char *name, int n, double *x)
{
  int i;
  for(i=0; i<n; i++){
    printf("%s[%d]=%+.15e\n",name,i,x[i]);
  }
}


void dqds_nakabayashi(int m, double *b, double *a, int *n_out, double *lambda,int debug)
{
  int i,k,done,nmax,splitflag,split,now,itr,kend,n1,n2,*mst=NULL,*width=NULL;
  double *q=NULL,*e=NULL,*q1=NULL,*q2=NULL,*e1=NULL,*W=NULL,Wmin,dmu,*d=NULL,*mu=NULL,*mu1=NULL,scale,x4,x5,x6,x7,eps2,wmin,qa,qb,ea,tmp,eps,s,t,m0,n;
  
  q=(double*)malloc(sizeof(double)*m);
  e=(double*)malloc(sizeof(double)*m-1);
  mst=(int*)malloc(sizeof(double)*m);
  mu=(double*)malloc(sizeof(double)*m);
  mu1=(double*)malloc(sizeof(double)*m);
  q1=(double*)malloc(sizeof(double)*m);
  e1=(double*)malloc(sizeof(double)*m-1);
  d=(double*)malloc(sizeof(double)*m);
  width=(int*)malloc(sizeof(double)*m);
  W=(double*)malloc(sizeof(double)*m);
  m0=m;
  eps=pow(2,-53);
  scale=pow(2,0);
  nmax=m*10;
  n=0;done=0;
  mu[0]=0;
  mst[0]=0;
  now=0;
  width[0]=m;
  itr=0;
  n1=0;
  n2=0;
  splitflag=1;
  //split_flag=0;
  eps2=eps*100*eps*100;
  
  
  for(int i=0;i<m;i++)
    {
      a[i]=a[i]*scale;
      b[i]=b[i]*scale;   
    }
  
  //LR分解　ズラシタ
  //aの要素単位のべき乗　
  for(i=0;i<m;i++)
    {
      a[i]=a[i]*a[i];
    }
  
  k=0;
  q[k]=b[k];
  //  print_vector("q",m,q);
  //  exit(0);
  k=k+1;
  while(k<m){
    e[k-1]=a[k-1]/q[k-1];
    q[k]=b[k]-e[k-1];
    k=k+1;
  }
  //print_vector("e",m-1,e);
  //print_vector("q",m,q);
  //exit(0);
  
  //////////////
  //終了フラグがオフかつ最大反復回数までループ　ズラシタ
  while(done==0 && n<=nmax){
    //printf("here2\n");
    split=1;
    while(split){if(splitflag==1){
	kend=mst[now]+1;
      }else{kend=MAX(mst[now]+1,m-2);}
      k=m-1;
      //1周目はe[98]=0.009183,q[99]=3.457618e-15,であってる
      
      while(k>=kend){
	if( fabs(e[k-1]) < fabs(q[k])*eps){
	  e[k-1]=0;
	  mst[now+1]=k;
	  width[now]=k-mst[now];
	  width[now+1]=m-k;
	  mu[now+1]=mu[now];
	  // printf("[%d~%d] -->[%d~%d]+[%d~%d]\n",mst[now],m-1,mst[now],mst[now+1]-1,mst[now+1],m-1);
	  now =now+1;
	  break;
	}
	k=k-1;
      }
      
      //１*１行列ズラシタ
      if(width[now]==1){
	lambda[m-1]=mu[now]+q[m-1];
	m=m-1;
	now = now-1;
	itr=0;
	if(m==0){
	  done=1;
	  split=0;
	  break;
	}
	continue;
      }
      
      //2*2行列
      //printf("now=%d \n",now);
      //printf("m=%d \n",m);
      //printf("q[%d]=%f \n",m-2,q[m-2]);
      //printf("q[%d]=%f \n",m-1,q[m-1]);
      //printf("e[%d]=%f \n",m-2,e[m-2]);
      //printf("here11a\n");fflush(stdout);
      else if(width[now]==2)
	{
	  //printf("here11\n");fflush(stdout);
	  //exit(0);
	  qa=q[m-2];
	  qb=q[m-1];
	  ea=e[m-2];
	  printf("qa=%f \n",qa);
	  //printf("qb=%f \n",qb);
	  //printf("ea=%f \n",ea);
	  
	  if(qa<qb)
	    {
	      tmp=qa;
	      qa=qb;
	      qb=tmp;
	    }
	  if(ea>eps2){
	    t=((qa-qb)+ea)/2;
	    //printf("t=%f \n",t);exit(0);
	    s=qb*(ea/t);
	    //printf("s=%f \n",s);exit(0);
	    if(s<=t)
	      {
		s=qb*ea/(t*(1+sqrt(1+s/t)));
	      }else{
	      s=qb*ea/(t+sqrt(t*(t+s)));
	    }
	    t=qa+(s+ea);
	    qb=qb*(qa/t);
	    qa=t;
	  }
	  lambda[m-1]=qa+mu[now];
	  //printf("lambda[%d]=%f \n",m-1,lambda[m-1]);
	  //lambda[1]=52.000000になった。MATLABではlambda[1]＝52.3970
	  lambda[m-2]=qb+mu[now];
	  //printf("lambda[%d]=%f \n",m-2,lambda[m-2]);
	  //lambda[0]=47.000000になった。MATLABではlambda[0]＝46.9423
	  m=m-2;
	  now=now-1;
	  itr=0;
	  if(m==0){
	    done=1;
	    split=0;
	    break;
	  }
	  continue;
	}
      split=0;break;
    }
    //exit(0);
    
    if(m==0){break;}
    //printf("here4\n");
    
    
    //シフト推定部分
    
    k=mst[now];
    W[k]=sqrt(q[k])-sqrt(e[k])/2;
    //printf("W[%d]=%f \n",k,W[k]); //１ループ目は問題なし。W[0]=7.179155
    k=k+1;
    while(k<m-1)
      {
	//printf("here12\n");fflush(stdout);
	W[k]=sqrt(q[k])-(sqrt(e[k-1])+sqrt(e[k]))/2;
	k=k+1;
      }
    //printf("k=%d \n",k);exit(0);
    
    W[k]=sqrt(q[k])-(sqrt(e[k-1]))/2;
    //printf("W[%d]=%f \n",k,W[k]);exit(0); //１ループ目は問題なし。W[99]=5.532719
    //printf("W[%d]=%f \n2",m-1,W[m-1]);exit(0);
    //１ループ目は問題なし。W[0]=7.179155、W[99]=5.53271
    Wmin=W[mst[now]];
    for (i=mst[now]; i<m; i++){ if (Wmin>W[i]){ Wmin=W[i]; } }
    if (Wmin<=0){ Wmin=0; }
    wmin=Wmin*Wmin; t=MAX(0,wmin); dmu=t;      
    
    //次のステップの計算
    k=mst[now];
    //printf("k=%d \n2",k);exit(0);
    d[k]=q[k]-dmu;
    //printf("d[%d]=%f \n2",k,d[k]);exit(0);
    //１ループ目は問題なし。d[0]=51.741492
    while(k<m){
      if(k<m-1){
	q1[k]=d[k]+e[k];
      }else if (k==m-1){
	q1[k]=d[k];
      }
      if(k<m-1){
	e1[k]=(q[k+1]/q1[k])*e[k];
	d[k+1]=(q[k+1]/q1[k])*d[k]-dmu;
      }
      k=k+1;
    }
    
    
    if (has_non_positive(m-mst[now],&q[mst[now]])||has_nan_or_inf(m-mst[now],&q1[mst[now]])||has_non_positive(m-mst[now]-1,&e[mst[now]])||has_nan_or_inf(m-mst[now]-1,&e1[mst[now]])){
      dmu=0;
      //次のステップ
      k=mst[now];
      d[k]=q[k]-dmu;
      while (k<m){
	if(k<m-1){
	  q1[k]=d[k]+e[k];
	}else if(k==m-1){
	  q1[k]=d[k];
	}
	if(k<m-1){
	  e1[k]=(q[k+1]/q1[k])*e[k];
	  d[k+1]=(q[k+1]/q1[k])*d[k]-dmu;
	}
	k=k+1;
      }
    }
    
    //printf("k=%d \n",k);
    //printf("q1[%d]=%f \n",k-1,q1[k-1]);
    //q1[99]=30.591137 
    //printf("e1[%d]=%f \n",k-2,e1[k-2]);
    //e1[98]=0.002914
    //printf("d[%d]=%f \n",k-1,d[k-1]);
    //d[99]=30.591137
    //exit(0);
    
    if(dmu==0){
      n1=n1+1;
    }else{
      n2=n2+1;
    }
    
    //正確にはmst(now)~mまで
    for(i=mst[now];i<m;i++){
      q[i]=q1[i];
    }
    for(i=mst[now];i<m-1;i++){
      e[i]=e1[i];
    }
    
    //printf("q[%d]=%f \n",m-1,q[m-1]);
    //printf("e[%d]=%f \n",m-2,e[m-2]);
    //exit(0);
    //q[99]=30.591137,e[98]=0.002914 
    
    //固有地の原点シフト量
    mu1[now]=mu[now]+dmu;
    mu[now]=mu1[now];
    //printf("mu1[%d]=%f \n",now,mu1[now]);
    //printf("mu[%d]=%f \n",now,mu[now]);
    //mu1[0]=0.549401,mu[0]=0.549401 
    //exit(0);
    n=n+1;
    itr=itr+1;
  }
  
  

    for(int i=0;i<=m;i++)
    {
      lambda[i]=lambda[i]/scale;   
    }
    
    printf("done=%d\n",done);
  //exit(0);
  printf("総反復回数: %f\n",n);
  //exit(0);
  double d1;
  d1=n/m0;
  printf("固有値一戸当たりの反復回数に平均: %f\n",d1);
  //exit(0);
  /*    
	for(i=0;i<n;i++){
	printf("固有値:%f \n",lambda[i]);
    }
  */
  

  //exit(0);



  //メモリ解放
  free(q);q=NULL;
  free(e);e=NULL;
  free(d);d=NULL;
  free(q1);q1=NULL;
  free(e1);e1=NULL;
  free(mst);mst=NULL;
  free(width);width=NULL;
  free(mu);mu=NULL;
  free(mu1);mu1=NULL; 

}



////////////////

int main(int argc, char *argv[])
{
  int i,j;
  FILE *fid=NULL;
  char type[3];   // ??? type="DST"
  int32_t n;      // 行列のサイズ
  double *A=NULL; // 行列の配列(n,n) -> n*n
  double *lambda=NULL; // 固有値の配列
  int n_out; // 固有値の計算できた個数
  double *b=NULL,*a=NULL; // 対角成分，．．
  int debug=1;
  clock_t start,end;
  double diff;

  // input 
  if((fid=fopen(FILENAME_INPUT,"r"))==NULL){
    printf("ファイルが見つかりません:FILENAME_INPUT¥n");
    exit(1);
  } // open  
  fread(type,sizeof(char),3,fid); // load
  fread(&n,sizeof(int32_t),1,fid); // load
  A=(double*)malloc(sizeof(double)*n*n); // allocate
  fread(A,sizeof(double),n*n,fid); // load
  fclose(fid);  

  //行列Aのチェック
/*
 for(j=0; j<n; j++){
    printf("\n");
    printf("j=%d\n",j);
    for(i=0; i<n; i++){
      printf("%+.3e\n",MAT(A,i,j,n));
    }
  }
*/

// 行列のサイズ
  //printf("n=%d\n",n);

    
  // eig
  lambda=(double*)malloc(sizeof(double)*n); // allocate
  b=(double*)malloc(sizeof(double)*n); // allocate
  a=(double*)malloc(sizeof(double)*(n-1)); // allocate
  for(i=0; i<n; i++){ b[i]=MAT(A,i,i,n); } // init
  for(i=0; i<n-1; i++){ a[i]=MAT(A,i,i+1,n); } // init
  for(i=0; i<n; i++){ lambda[i]=0; }        // initialize
  

//  print_vector("b",n,b);
//  print_vector("a",n-1,a);
//  print_vector("lambda",n,lambda);

  // dqds
  start=clock();
  dqds_nakabayashi(n,b,a,&n_out,lambda,debug);
  end=clock();
  //  printf("n_out=%d\n",n_out);
  print_vector("lambda",n,lambda);
  diff=(double)(end-start)/CLOCKS_PER_SEC;
  printf("%.15e [s] \n",diff);
  
  //exit(0);
  // output
  if((fid=fopen(FILENAME_OUTPUT,"w"))==NULL){
    printf("ファイルが見つかりません:FILENAME_OUTPUT¥n");
    exit(1);
  }// open
  //exit(0);
  fwrite(type,sizeof(char),3,fid); // save
  fwrite(&n,sizeof(int32_t),1,fid); // save
  fwrite(lambda,sizeof(double),n,fid); // save
  fwrite(&diff,sizeof(double),1,fid);
  fclose(fid);

  // done
  free(A); A=NULL;
  free(lambda); lambda=NULL;
  free(b); b=NULL;
  free(a); a=NULL;
  return 0;
}
