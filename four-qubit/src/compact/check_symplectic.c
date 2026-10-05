/* Exhaustive binary implications: common ordered chart or second Pauli return. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
typedef unsigned char U;
typedef struct {U g[6]; unsigned n; unsigned *ids;} St;
typedef struct {uint64_t mask;U bs[3];} Edge;
typedef struct {unsigned nt,ns,ne,nr,nb;U (*t)[6];St *s;Edge *e;U (*r)[2],(*b)[2];uint64_t pre,common,ret,bal,filtered,bad;U first[6];} Case;
static U F[6],om[64][64],basis[135][3];static uint64_t lm[135],total,uni;static Case *cases;static unsigned nc;
static void die(const char*s){fprintf(stderr,"%s\n",s);exit(2);}
static void rd(FILE*f,void*p,size_t n){if(fread(p,1,n,f)!=n)die("truncated");}
static void six(const U*a,unsigned n){for(unsigned i=0;i<n;i++)if(a[i]>=64)die("not six bit");}
static unsigned ui(FILE*f){unsigned x;rd(f,&x,4);return x;}
static U app(const U*a,U v){U r=0;for(int j=0;j<6;j++)if(v&(1<<j))r^=a[j];return r;}
static void mul(U*r,const U*a,const U*b){for(int j=0;j<6;j++)r[j]=app(a,b[j]);}
static int unip(const U*m){U a[6],b[6];for(int j=0;j<6;j++)a[j]=m[j]^(1<<j);for(int k=0;k<3;k++){mul(b,a,a);memcpy(a,b,6);}for(int j=0;j<6;j++)if(a[j])return 0;return 1;}
static int pairs(U (*p)[2],unsigned n){for(unsigned j=0;j<n;j++)if(app(F,p[j][1])==p[j][0])return 1;return 0;}
static void visit(void){
 total++;if(!unip(F))return;uni++;U a[6];
 for(unsigned k=0;k<nc;k++){
  Case*c=&cases[k];int ok=1;
  for(unsigned j=1;j<c->nt;j++){mul(a,F,c->t[j]);if(!unip(a)){ok=0;break;}}
  if(!ok)continue;c->pre++;
  for(unsigned j=0;j<c->ne;j++){
   Edge*e=&c->e[j];
   if((e->mask>>app(F,e->bs[0])&1)&&(e->mask>>app(F,e->bs[1])&1)&&(e->mask>>app(F,e->bs[2])&1)){ok=2;break;}
  }
  if(ok==2){c->common++;continue;}
  if(pairs(c->r,c->nr)){c->ret++;continue;}
  if(!c->bad)memcpy(c->first,F,6);c->bad++;
 }
}
static void enumerate(int d){
 if(d==3){visit();return;}
 for(int a=1;a<64;a++){
  int ok=1;for(int j=0;j<d;j++)if(om[a][F[j]]||om[a][F[j+3]])ok=0;
  if(!ok)continue;F[d]=a;
  for(int b=1;b<64;b++){
   if(!om[a][b])continue;ok=1;for(int j=0;j<d;j++)if(om[b][F[j]]||om[b][F[j+3]])ok=0;
   if(ok){F[d+3]=b;enumerate(d+1);}
  }
 }
}
int main(int argc,char**argv){
 if(argc!=2)die("need data path");
 for(int a=0;a<64;a++)for(int b=0;b<64;b++){int v=((a&7)&(b>>3))^((b&7)&(a>>3));om[a][b]=(v^(v>>1)^(v>>2))&1;}
 FILE*f=fopen(argv[1],"rb");if(!f)die("no data");unsigned nl=ui(f);if(nl!=135)die("lag count");
 for(unsigned l=0;l<135;l++){
  rd(f,&lm[l],8);rd(f,basis[l],3);six(basis[l],3);uint64_t mask=0;
  for(int t=0;t<8;t++){int v=0;for(int j=0;j<3;j++)if(t>>j&1)v^=basis[l][j];mask|=(uint64_t)1<<v;}
  if(mask!=lm[l]||__builtin_popcountll(mask)!=8)die("lag span");
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(om[basis[l][i]][basis[l][j]])die("lag isotropy");
  for(unsigned j=0;j<l;j++)if(lm[j]==lm[l])die("duplicate lag");
 }
 nc=ui(f);if(!nc||nc>100000)die("case count");cases=calloc(nc,sizeof(Case));
 for(unsigned k=0;k<nc;k++){
  Case*c=&cases[k];c->nt=ui(f);c->ns=ui(f);c->ne=ui(f);c->nr=ui(f);c->nb=ui(f);
  if(!c->nt||c->nt>4096||c->ns>64||c->ne>135||c->nr>4096||c->nb>4096)die("case size");
  c->t=calloc(c->nt,6);rd(f,c->t,6*c->nt);six((U*)c->t,6*c->nt);
  for(int j=0;j<6;j++)if(c->t[0][j]!=(1<<j))die("identity prefilter");
  c->s=calloc(c->ns,sizeof(St));
  for(unsigned j=0;j<c->ns;j++){St*s=&c->s[j];rd(f,s->g,6);six(s->g,6);s->n=ui(f);if(s->n>135)die("station size");s->ids=calloc(s->n,4);for(unsigned h=0;h<s->n;h++){s->ids[h]=ui(f);if(s->ids[h]>=135)die("station index");}}
  c->e=calloc(c->ne,sizeof(Edge));for(unsigned j=0;j<c->ne;j++){rd(f,&c->e[j].mask,8);rd(f,c->e[j].bs,3);six(c->e[j].bs,3);}
  c->r=calloc(c->nr,2);rd(f,c->r,2*c->nr);six((U*)c->r,2*c->nr);c->b=calloc(c->nb,2);rd(f,c->b,2*c->nb);six((U*)c->b,2*c->nb);
 }
 if(fgetc(f)!=EOF)die("trailing data");fclose(f);clock_t start=clock();enumerate(0);
 printf("{\"symplectic_matrices\":%llu,\"unipotent\":%llu,\"seconds\":%.6f,\"cases\":[",(unsigned long long)total,(unsigned long long)uni,(double)(clock()-start)/CLOCKS_PER_SEC);
 uint64_t bad=0;
 for(unsigned k=0;k<nc;k++){Case*c=&cases[k];bad+=c->bad;
 printf("%s{\"id\":%u,\"unipotence_pass\":%llu,\"common\":%llu,\"returned\":%llu,\"balanced\":%llu,\"stationarity_excluded\":%llu,\"bad\":%llu,\"first\":[%u,%u,%u,%u,%u,%u]}",k?",":"",k,(unsigned long long)c->pre,(unsigned long long)c->common,(unsigned long long)c->ret,(unsigned long long)c->bal,(unsigned long long)c->filtered,(unsigned long long)c->bad,c->first[0],c->first[1],c->first[2],c->first[3],c->first[4],c->first[5]);}
 printf("],\"bad\":%llu}\n",(unsigned long long)bad);return bad||total!=1451520||uni!=262144;
}
