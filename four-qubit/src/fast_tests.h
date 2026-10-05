/* Independent brute-force counts against exact Z/4 elimination and family cover. */
static uint32_t rng_state=20260915;
static uint32_t frand(void){rng_state^=rng_state<<13;rng_state^=rng_state>>17;rng_state^=rng_state<<5;return rng_state;}
static uint64_t brute_count(const fc_system *s){
    uint64_t count=0,all=UINT64_C(1)<<(2*s->n);
    for(uint64_t z=0;z<all;z++){
        int good=1;
        for(int i=0;i<s->m && good;i++){
            unsigned sum=0;for(int j=0;j<s->n;j++)sum+=s->a[i][j]*((z>>(2*j))&3);
            if((sum&3)!=s->a[i][16])good=0;
        }
        count+=good;
    }
    return count;
}
static int dense_word_meets(const int *pi,const int *f,int lam) {
    for(int k=0;k<3;k++){
        int ray=fast_pulled[k],p=ray&15,sign=-1,ok=1;
        for(int x=0;x<N;x++){
            if(pi[x^p]!=(pi[x]^p)){ok=0;break;}
            int left=(pw(ray,x)+lam+2*f[x^p])&7;
            int right=(lam+2*f[x]+pw(ray,pi[x]))&7;
            int ratio=(left-right)&7;
            if(ratio!=0 && ratio!=4){ok=0;break;}
            if(sign<0)sign=ratio;else if(sign!=ratio){ok=0;break;}
        }
        if(ok)return 1;
    }
    return 0;
}

static int fast_tests(void){
    for(int it=0;it<12000;it++){
        fc_system s;memset(&s,0,sizeof s);s.n=frand()%7;s.m=frand()%19;
        for(int i=0;i<s.m;i++){for(int j=0;j<s.n;j++)s.a[i][j]=frand()%4;s.a[i][16]=frand()%4;}
        uint64_t x=fc_count(s),y=brute_count(&s);
        if(x!=y){printf("COUNT_FAIL %d %llu %llu\n",it,(unsigned long long)x,(unsigned long long)y);return 1;}
    }
    /* Identity point map with arbitrary independent diagonal signs is Hermitian;
       test affine subfamilies and arbitrary planes, including variable-ray covers. */
    int pi[N];for(int x=0;x<N;x++)pi[x]=x;
    for(int it=0;it<2000;it++){
        ff_n=1+frand()%5;
        for(int x=0;x<N;x++){
            ff_base[x]=2*(frand()%2);
            for(int j=0;j<N;j++)ff_coef[x][j]=j<ff_n?2*(frand()%2):0;
        }
        fast_pulled[0]=1+frand()%255;
        do{fast_pulled[1]=1+frand()%255;}while(fast_pulled[0]==fast_pulled[1]);
        fast_pulled[2]=fast_pulled[0]^fast_pulled[1];
        uint64_t all=UINT64_C(1)<<(2*ff_n);int every=1;
        for(uint64_t z=0;z<all;z++){
            int f[N];for(int x=0;x<N;x++){
                int v=ff_base[x];for(int j=0;j<ff_n;j++)v+=ff_coef[x][j]*((z>>(2*j))&3);f[x]=v&3;
            }
            if(fast_meets(pi,f)!=dense_word_meets(pi,f,0)){printf("RAY_FORMULA_FAIL\n");return 1;}
            if(!dense_word_meets(pi,f,0)){every=0;break;}
        }
        int got=fc_cover(pi);if(got!=every){printf("FAMILY_FAIL %d\n",it);return 1;}
    }
    for(int it=0;it<160;it++){
        int order[N];for(int x=0;x<N;x++)order[x]=x;
        for(int x=15;x>0;x--){int y=frand()%(x+1),v=order[x];order[x]=order[y];order[y]=v;}
        int lam=it&1,kappa=lam?3:0;ff_n=1+frand()%5;
        memset(ff_coef,0,sizeof ff_coef);
        for(int k=0;k<8;k++){
            int x=order[2*k],y=order[2*k+1],j=frand()%ff_n;
            pi[x]=y;pi[y]=x;ff_base[x]=frand()%4;ff_base[y]=(kappa-ff_base[x])&3;
            ff_coef[x][j]=1;ff_coef[y][j]=3;
        }
        fast_pulled[0]=1+frand()%255;
        do{fast_pulled[1]=1+frand()%255;}while(fast_pulled[0]==fast_pulled[1]);
        fast_pulled[2]=fast_pulled[0]^fast_pulled[1];
        int every=1;uint64_t all=UINT64_C(1)<<(2*ff_n);
        for(uint64_t z=0;z<all;z++){
            int f[N];for(int x=0;x<N;x++){
                int v=ff_base[x];for(int j=0;j<ff_n;j++)v+=ff_coef[x][j]*((z>>(2*j))&3);f[x]=v&3;
            }
            int actual=dense_word_meets(pi,f,lam);
            if(fast_meets(pi,f)!=actual){printf("ODD_RAY_FAIL\n");return 1;}
            if(!actual)every=0;
        }
        if(fc_cover(pi)!=every){printf("ODD_FAMILY_FAIL\n");return 1;}
    }
    printf("FAST_TESTS_PASS modular_systems=12000 affine_families=2160\n");return 0;
}

