/* Exact phase-family tests used by the census kernel; the regression test fasttests compares them
   with brute force (fast_tests.h). */
static int fast_enabled=0, fast_counts=1, fast_first_chart=0, fast_last_chart=NLAG-1;
static int fast_plane[3], fast_pulled[3], fast_id, fast_cert_mode[1100], fast_cert_plane[1100][3];
static long long fast_early=0, fast_common=0, fast_four=0, fast_phase_skip=0, fast_quick=0;
static int fast_family_ready=0;
typedef struct { int n,m; unsigned char a[64][17]; } fc_system;
static uint64_t fc_count(fc_system s) {
    int r=0, n=s.n;
    for (;;) {
        int pr=-1, pc=-1;
        for(int i=r;i<s.m && pr<0;i++) for(int j=r;j<n;j++)
            if(s.a[i][j]&1){pr=i;pc=j;break;}
        if(pr<0) break;
        for(int j=0;j<=16;j++){unsigned char t=s.a[r][j];s.a[r][j]=s.a[pr][j];s.a[pr][j]=t;}
        for(int i=0;i<s.m;i++){unsigned char t=s.a[i][r];s.a[i][r]=s.a[i][pc];s.a[i][pc]=t;}
        if(s.a[r][r]==3) {
            for(int j=r;j<n;j++) s.a[r][j]=(3*s.a[r][j])&3;
            s.a[r][16]=(3*s.a[r][16])&3;
        }
        for(int i=r+1;i<s.m;i++){
            int v=s.a[i][r];
            for(int j=r;j<n;j++) s.a[i][j]=(s.a[i][j]-v*s.a[r][j])&3;
            s.a[i][16]=(s.a[i][16]-v*s.a[r][16])&3;
        }
        r++;
    }
    int rank1=r;
    for(int i=r;i<s.m;i++){
        if(s.a[i][16]&1) return 0;
        for(int j=r;j<n;j++) s.a[i][j]>>=1;
        s.a[i][16]>>=1;
    }
    int rr=r;
    for(int col=r;col<n;col++){
        int pivot=rr;while(pivot<s.m && !s.a[pivot][col])pivot++;
        if(pivot==s.m)continue;
        for(int j=0;j<=16;j++){unsigned char t=s.a[rr][j];s.a[rr][j]=s.a[pivot][j];s.a[pivot][j]=t;}
        for(int i=rr+1;i<s.m;i++)if(s.a[i][col]){
            for(int j=col;j<n;j++)s.a[i][j]^=s.a[rr][j];
            s.a[i][16]^=s.a[rr][16];
        }
        rr++;
    }
    for(int i=rr;i<s.m;i++)if(s.a[i][16])return 0;
    return UINT64_C(1)<<(2*(n-rank1)-(rr-rank1));
}
static int ff_n,ff_base[N],ff_coef[N][N];
static void fc_prepare(const int *alpha,const int *sgn,const int *part,const int *twice,int npart) {
    ff_n=npart; memset(ff_coef,0,sizeof ff_coef);
    for(int x=0;x<N;x++){
        int c=comp_of[x],p=part[c],offset=twice[p]<0?0:twice[p]/2,step=twice[p]<0?1:2;
        ff_base[x]=(alpha[c]+delta_of[x]+sgn[c]*offset)&3;
        ff_coef[x][p]=(sgn[c]*step)&3;
    }
    fast_family_ready=1;
}
static int fc_add_ray(fc_system *s,const int *pi,int ray) {
    int p=ray&15,q=ray>>4;
    for(int x=0;x<N;x++)if(pi[x^p]!=(pi[x]^p))return 0;
    int k0=(ff_base[p]-ff_base[0]+2*(popc(q&pi[0])&1))&3;
    for(int x=1;x<N;x++){
        int row=s->m++;
        int kx=(ff_base[x^p]-ff_base[x]+2*(popc(q&(pi[x]^x))&1))&3;
        s->a[row][16]=(k0-kx)&3;
        for(int j=0;j<ff_n;j++)
            s->a[row][j]=(ff_coef[x^p][j]-ff_coef[x][j]-ff_coef[p][j]+ff_coef[0][j])&3;
    }
    int row=s->m++; s->a[row][16]=(-2*k0)&3;
    for(int j=0;j<ff_n;j++)s->a[row][j]=(2*(ff_coef[p][j]-ff_coef[0][j]))&3;
    return 1;
}
static int fc_universal(const fc_system *s) {
    for(int i=0;i<s->m;i++){
        if(s->a[i][16])return 0;
        for(int j=0;j<s->n;j++)if(s->a[i][j])return 0;
    }
    return 1;
}
static int fc_cover(const int *pi) {
    fc_system systems[3]; uint64_t counts[3]={0,0,0};
    uint64_t all=UINT64_C(1)<<(2*ff_n);
    int valid[3];
    for(int k=0;k<3;k++){
        memset(&systems[k],0,sizeof(fc_system));systems[k].n=ff_n;
        valid[k]=fc_add_ray(&systems[k],pi,fast_pulled[k]);
        if(valid[k] && fc_universal(&systems[k])){fast_common++;return 1;}
    }
    for(int k=0;k<3;k++)if(valid[k])counts[k]=fc_count(systems[k]);
    uint64_t both=0;
    if(valid[0] && valid[1]){
        fc_system joint=systems[0];
        for(int i=0;i<systems[1].m;i++)memcpy(joint.a[joint.m++],systems[1].a[i],17);
        both=fc_count(joint);
    }
    if(counts[0]+counts[1]+counts[2]==all+2*both){fast_four++;return 1;}
    return 0;
}
static int fast_generic(const int *pi) {
    for(int k=0;k<3;k++){
        int r=fast_pulled[k]; if(r&15)continue;
        int q=r>>4,v=popc(q&pi[0])&1,ok=1;
        for(int x=1;x<N;x++)if((popc(q&(pi[x]^x))&1)!=v){ok=0;break;}
        if(ok)return 1;
    }
    return 0;
}
static int fast_meets(const int *pi,const int *f) {
    for(int z=0;z<3;z++){
        int r=fast_pulled[z],p=r&15,q=r>>4,ok=1;
        int k0=(f[p]-f[0]+2*(popc(q&pi[0])&1))&3;
        if(k0&1)continue;
        for(int x=0;x<N;x++)
            if(pi[x^p]!=(pi[x]^p) || ((f[x^p]-f[x]+2*(popc(q&(pi[x]^x))&1))&3)!=k0){ok=0;break;}
        if(ok)return 1;
    }
    return 0;
}
static void fast_load_certificate(void) {
    FILE *f=fopen("lock-n4.txt","r");if(!f){perror("lock-n4.txt");exit(1);}
    char line[512];int count=0;
    while(fgets(line,sizeof line,f)){
        int id,level,sector,n,a,b,c;char mode[32];
        if(sscanf(line,"%d level %d sector %d N %d mode %31s",&id,&level,&sector,&n,mode)!=5 ||
           id<0 || id>=1020 || fast_cert_mode[id]){fprintf(stderr,"invalid certificate row\n");exit(1);}
        if(!strcmp(mode,"largeN"))fast_cert_mode[id]=3;
        else {
            char *p=strstr(line,"plane (");
            if(!p || sscanf(p,"plane (%d,%d,%d)",&a,&b,&c)!=3 || a<=0 || a>=256 || b<=0 || b>=256 ||
                a==b || c!=(a^b)){fprintf(stderr,"invalid plane\n");exit(1);}
            fast_cert_mode[id]=!strcmp(mode,"pair")?1:!strcmp(mode,"exact")?2:0;
            if(!fast_cert_mode[id])exit(1);
            fast_cert_plane[id][0]=a;fast_cert_plane[id][1]=b;fast_cert_plane[id][2]=c;
        }
        count++;
    }
    fclose(f);if(count!=1020){fprintf(stderr,"incomplete proposed certificate\n");exit(1);}
}

static long long fast_rigid_skip=0;
static int fast_subset=0, fast_chart_selected[NLAG];
static int fast_color[N][N],fast_entry_phase[N][N];
static void fast_colorize(z8 A[N][N]){
    z8 representatives[256];int n=0;
    for(int x=0;x<N;x++)for(int y=0;y<N;y++){
        if(ziszero(A[x][y])){fast_color[x][y]=-1;fast_entry_phase[x][y]=0;continue;}
        int found=-1,phase=0;
        for(int j=0;j<n && found<0;j++)for(int r=0;r<4;r++)
            if(zeq(A[x][y],zmul(zunit(2*r),representatives[j]))){found=j;phase=r;break;}
        if(found<0){found=n;representatives[n++]=A[x][y];}
        fast_color[x][y]=found;fast_entry_phase[x][y]=phase;
    }
}

