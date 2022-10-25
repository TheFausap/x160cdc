//
//  main.c
//  x160cdc
//
//  Created by Fausto Saporito on 20/10/2022.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ovr;
int is_hlt;
int BFR_busy = 0;
int dvcstatus = 0;

int* A;
int* Z;
int* P;
int* G;
int* S;
int* F;
int* E;
int* Ap;

int B;
int I;
int R;
int D;

int* PSR;
int* BER;
int* BFR;
int* BXR;

int* DW;

int* sbc;

int*** mem;

#define RS 12

int on1[] = {0,0,0,0,0,0,0,0,0,0,0,1}; // +1
int zr0[] = {0,0,0,0,0,0,0,0,0,0,0,0}; // +0
int zr1[] = {1,1,1,1,1,1,1,1,1,1,1,1}; // -0

int sw1;
int sw2;
int sw4;

FILE* tape;
FILE* prt;
FILE* ptape;

FILE* dvc;

void _set(int* d, int s)
{
    int i = 11;
    
    while(s != 0) {
        d[i] = s%2;
        s /= 2;
        i--;
    }
}

int hs(int a, int b, int* bo)
{
    (*bo) = (!a) & b;
    return a ^ b;
}

int fs(int a, int b, int bi, int* bo)
{
    int d = 0;
    int bo1 = 0;
    int bo2 = 0;
    
    d = hs(hs(a,b,&bo1),bi,&bo2);
    (*bo) = bo1 | bo2;
    
    return d;
}

void cp(int* a)
{
    for(int i=0;i<RS;i++) {
        a[i] = !a[i];
    }
}


void lp(int* a,int *b)
{
    for(int i=0;i<RS;i++) {
        a[i] = a[i] & b[i];
    }
}

// copy lower (E part)
void cpl(int* d, int* s)
{
    for(int i=0;i<6;i++) {
        d[6+i] = s[i];
    }
}

// copy upper (F part)
void cpu(int*d, int* s)
{
    for(int i=0;i<6;i++) {
        d[i]=s[i];
    }
}

// read one word from the device d
int* rdd(FILE* d)
{
    int* t;
    t=calloc(12,sizeof(int));
    for(int i=0;i<RS;i++) {
        t[i] = fgetc(d)-'0';
    }
    
    return t;
}

// write one word to the device d
void wrd(FILE* d, int* s)
{
    for(int i=0;i<RS;i++) {
        fputc(s[i]+'0',d);
    }
}

void adder(int* a, int* b, int op){
    int bi = 0;
    int bc[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for (int i=0;i<RS;i++) {
        bc[i] = b[i];
    }
    
    if (op)
        cp(bc);
    for(int i=RS-1;i>=0;i--) {
        a[i] = fs(a[i],bc[i],bi,&ovr);
        bi = ovr;
    }
    if (ovr == 1) {
        for(int i=RS-1;i>=0;i--) {
            a[i] = fs(a[i],zr0[i],bi,&ovr);
            bi = ovr;
        }
    }
}

void inc(int* a){
    int bi = 0;
    int one[] = {0,0,0,0,0,0,0,0,0,0,0,1};
    
    cp(one);
    for(int i=RS-1;i>=0;i--) {
        a[i] = fs(a[i],one[i],bi,&ovr);
        bi = ovr;
    }
    if (ovr == 1) {
        for(int i=RS-1;i>=0;i--) {
            a[i] = fs(a[i],zr0[0],bi,&ovr);
            bi = ovr;
        }
    }
}

int ne(int* a, int* b)
{
    for (int i=0;i<RS;i++) {
        if (a[i] != b[i]) {
            return 1;
        }
    }
    return 0;
}

void init(void)
{
    FILE* pa;
    int t;
    
    A = calloc(12,sizeof(int));
    P = calloc(12,sizeof(int));
    BER = calloc(12,sizeof(int));
    BXR = calloc(12,sizeof(int));
    
    DW = calloc(12,sizeof(int));
    
    BFR = calloc(12,sizeof(int));
    S = calloc(12,sizeof(int));
    F = calloc(6,sizeof(int));
    E = calloc(6,sizeof(int));
    G = calloc(12,sizeof(int));
    Ap = calloc(12,sizeof(int));
    PSR = calloc(8,sizeof(int));
    Z = calloc(12,sizeof(int));
    
    sbc = calloc(4,sizeof(int));
    
    mem = calloc(8,sizeof(int**));
    for (int i=0;i<8;i++) {
        mem[i] = calloc(4096,sizeof(int*));
        for (int j=0;j<4096;j++) {
            mem[i][j] = calloc(12,sizeof(int));
        }
    }
    
    B=0;I=0;R=0;D=0;
    sw1=0;sw2=0;sw4=0;
    
    pa = fopen("config.ini","r");
    
    if (pa == NULL) {
        printf("Proceeding without specific setup from the panel\n");
    } else {
        printf("Reading setup from the panel\n");
        // config file format
        // AAAA B I R D PPPP SW1 SW2 SW4
        // reading A
        t = (fgetc(pa)-'0')*512 + (fgetc(pa)-'0')*64 + (fgetc(pa)-'0')*8 +
            (fgetc(pa)-'0');
        _set(A,t);
        // reading B
        fgetc(pa);
        B = fgetc(pa)-'0';
        //reading I
        fgetc(pa);
        I = fgetc(pa)-'0';
        // reading R
        fgetc(pa);
        R = fgetc(pa)-'0';
        // reading D
        fgetc(pa);
        D = fgetc(pa)-'0';
        // reading P
        fgetc(pa);
        t = (fgetc(pa)-'0')*512 + (fgetc(pa)-'0')*64 + (fgetc(pa)-'0')*8 +
            (fgetc(pa)-'0');
        _set(P,t);
        // reading SW1
        fgetc(pa);
        sw1 = fgetc(pa)-'0';
        // reading SW2
        fgetc(pa);
        sw2 = fgetc(pa)-'0';
        // reading SW4
        fgetc(pa);
        sw4 = fgetc(pa)-'0';
    }
}

int to_num(int* b, int l)
{
    int r = 0;
    
    for(int i=l-1,j=0;i>=0;i--,j++) {
        r += b[i] * pow(2,j);
    }
    
    return r;
}

void set(int* d, int* s)
{
    for(int i=0;i<RS;i++)
        d[i] = s[i];
}

void from_mem(int bl, int* d, int* s, int l)
{
    for(int i=0;i<RS;i++) {
        d[i] = mem[bl][to_num(s,l)][i];
    }
}

void to_mem(int bl, int* d, int l, int* s)
{
    for(int i=0;i<RS;i++) {
        mem[bl][to_num(d,l)][i]=s[i];
    }
}

void readFE(void)
{
    int FE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    from_mem(R,FE,P,12);
    
    F[0]=FE[0];F[1]=FE[1];F[2]=FE[2];
    F[3]=FE[3];F[4]=FE[4];F[5]=FE[5];
    
    E[0]=FE[6];E[1]=FE[7];E[2]=FE[8];
    E[3]=FE[9];E[4]=FE[10];E[5]=FE[11];
}

void readG(void)
{
    inc(P);
    from_mem(R,G,P,12);
}

void dump(void);

void op_err(void)
{
    is_hlt = 1;
    dump();
    printf("\nERR\n");
}

void op_hlt(int t)
{
    is_hlt = 1;
    dump();
    if (t == 0) {
        printf("\n7700\n");
    } else {
        printf("\n7777\n");
    }
}

void op_nop(void)
{
    // nothing, nic!
}

void op_bls(void)
{
    if (BFR_busy) {
        set(P,G);
    } else {
        BFR_busy = 1;
        set(BFR,A);
        while (ne(BER,BXR)) {
            to_mem(B,BER,12,BFR);
            inc(BER);
        }
    }
    BFR_busy = 0;
}

void op_pta(void)
{
    set(A,P);
}

void op_ate(void)
{
    if (BFR_busy) {
        set(P,G);
    } else {
        set(BER,A);
    }
}

void op_atx(void)
{
    if (BFR_busy) {
        set(P,G);
    } else {
        set(BXR,A);
    }
}

void op_eta(void)
{
    set(A,BER);
}

void op_cta(void)
{
    int _B=B;int _D=D;int _I=I;int _R=R;
    
    A[11]=_B%2;_B/=2;A[10]=_B%2;_B/=2;A[9]=_B%2;
    A[8]=_D%2;_D/=2;A[7]=_D%2;_D/=2;A[6]=_D%2;
    A[5]=_I%2;_I/=2;A[4]=_I%2;_I/=2;A[3]=_I%2;
    A[2]=_R%2;_R/=2;A[1]=_R%2;_R/=2;A[0]=_R%2;
}

void op_stp(void)
{
    for(int i=0;i<11;i++) {
        mem[D][to_num(E,6)][i] = P[i];
    }
}

void op_ldn(void)
{
    for(int i=0;i<6;i++) {
        A[i] = 0;
    }
    for(int i=6;i<RS;i++) {
        A[i] = E[i-6];
    }
}

void op_ldd(void)
{
    for(int i=0;i<RS;i++) {
        A[i] = mem[D][to_num(E,6)][i];
    }
}

void op_ldi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    for(int i=0;i<RS;i++) {
        A[i] = mem[I][ad][i];
    }
}

void op_ldf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        A[i] = mem[R][ad][i];
    }
}

void op_ldb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        A[i] = mem[R][ad][i];
    }
}

void op_lds(void)
{
    for(int i=0;i<RS;i++) {
        A[i] = mem[0][07777][i];
    }
}

void op_ldc(void)
{
    set(A,G);
}

void op_ldm(void)
{
    set(A,mem[I][to_num(G,12)]);
}

void op_lcn(void)
{
    for(int i=0;i<6;i++) {
        A[i] = 1;
    }
    for(int i=6;i<RS;i++) {
        A[i] = !E[i-6];
    }
}

void op_lcd(void)
{
    for(int i=0;i<RS;i++) {
        A[i] = !mem[D][to_num(E,6)][i];
    }
}

void op_lci(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    for(int i=0;i<RS;i++) {
        A[i] = !mem[I][ad][i];
    }
}

void op_lcf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        A[i] = !mem[R][ad][i];
    }
}

void op_lcb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        A[i] = !mem[R][ad][i];
    }
}

void op_lcs(void)
{
    for(int i=0;i<RS;i++) {
        A[i] = !mem[0][07777][i];
    }
}

void op_lcc(void)
{
    set(A,G);
    cp(A);
}

void op_lcm(void)
{
    set(A,mem[I][to_num(G,12)]);
    cp(A);
}

void op_ste(void)
{
    to_mem(D,E,6,BER);
    set(BER,A);
}

void op_std(void)
{
    for(int i=0;i<RS;i++) {
        mem[D][to_num(E,6)][i] = A[i];
    }
}

void op_sti(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    for(int i=0;i<RS;i++) {
        mem[I][ad][i] = A[i];
    }
}

void op_stf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        mem[R][ad][i] = A[i];
    }
}

void op_stb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        mem[R][ad][i] = A[i];
    }
}

void op_sts(void)
{
    for(int i=0;i<RS;i++) {
        mem[0][07777][i] = A[i];
    }
}

void op_stc(void)
{
    set(G,A);
}

void op_stm(void)
{
    set(mem[I][to_num(G,12)],A);
}

void op_hwi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    for(int i=6;i<RS;i++) {
        mem[I][ad][i] = A[i];
    }
}

void op_ina(void)
{
    if (dvcstatus) {
        set(A,DW);
    } else {
        for(int i=0;i<RS;i++) {
            A[i] = fgetc(dvc)-'0';
        }
    }
}

void op_ota(void)
{
    for(int i=0;i<RS;i++) {
        fputc(A[i]+'0',dvc);
    }
}

void op_mut(void)
{
    for(int i=0;i<10;i++) {
        adder(A,A,1);
    }
}

void op_muh(void)
{
    for(int i=0;i<100;i++) {
        adder(A,A,1);
    }
}

void op_adn(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for(int i=6;i<RS;i++) {
        t[i] = E[i-6];
    }
    adder(A,t,1);
}

void op_add(void)
{
    adder(A,mem[D][to_num(E,6)],1);
}

void op_adi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    adder(A,mem[I][ad],1);
}

void op_adf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    adder(A,mem[R][ad],1);
}

void op_adb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    adder(A,mem[R][ad],1);
}

void op_ads(void)
{
    adder(A,mem[0][07777],1);
}

void op_adc(void)
{
    adder(A,G,1);
}

void op_adm(void)
{
    adder(A,mem[I][to_num(G,12)],1);
}

void op_sbn(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for(int i=6;i<RS;i++) {
        t[i] = E[i-6];
    }
    adder(A,t,0);
}

void op_sbd(void)
{
    adder(A,mem[D][to_num(E,6)],0);
}

void op_sbi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    adder(A,mem[I][ad],0);
}

void op_sbf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    adder(A,mem[R][ad],0);
}

void op_sbb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    adder(A,mem[R][ad],0);
}

void op_sbs(void)
{
    adder(A,mem[0][07777],0);
}

void op_sbc(void)
{
    adder(A,G,0);
}

void op_sbm(void)
{
    adder(A,mem[I][to_num(G,12)],0);
}

void op_rad(void)
{
    adder(A,mem[D][to_num(E,6)],1);
    to_mem(D,E,6,A);
}

void op_rai(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    adder(A,mem[I][ad],1);
    to_mem(I,E,6,A);
}

void op_raf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    adder(A,mem[R][ad],1);
    set(mem[R][ad],A);
}

void op_rab(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    adder(A,mem[R][ad],1);
    set(mem[R][ad],A);
}

void op_ras(void)
{
    adder(A,mem[0][07777],1);
    set(mem[0][07777],A);
}

void op_rac(void)
{
    adder(A,G,1);
    set(G,A);
}

void op_ram(void)
{
    adder(A,mem[I][to_num(G,12)],1);
    to_mem(I,G,12,A);
}

void op_aod(void)
{
    set(A,mem[D][to_num(E,6)]);
    inc(A);
    to_mem(D,E,6,A);
}

void op_aoi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    set(A,mem[I][ad]);
    inc(A);
    to_mem(I,E,6,A);
}

void op_aof(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    set(A,mem[R][ad]);
    inc(A);
    set(mem[R][ad],A);
}

void op_aob(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    set(A,mem[R][ad]);
    inc(A);
    set(mem[R][ad],A);
}

void op_aos(void)
{
    set(A,mem[0][07777]);
    inc(A);
    set(mem[0][07777],A);
}

void op_aoc(void)
{
    set(A,G);
    inc(A);
    set(G,A);
}

void op_aom(void)
{
    set(A,mem[I][to_num(G,12)]);
    inc(A);
    to_mem(I,G,12,A);
}

void op_ls1(void)
{
    int t = A[0];
    
    A[0]=A[1];A[1]=A[2];
    A[2]=A[3];A[3]=A[4];
    A[4]=A[5];A[5]=A[6];
    A[6]=A[7];A[7]=A[8];
    A[8]=A[9];A[9]=A[10];
    A[10]=A[11];A[11]=t;
}

void op_ls2(void)
{
    op_ls1();
    op_ls1();
}

void op_ls3(void)
{
    op_ls1();
    op_ls1();
    op_ls1();
}

void op_ls6(void)
{
    op_ls1();
    op_ls1();
    op_ls1();
    op_ls1();
    op_ls1();
    op_ls1();
}

void op_rs1(void)
{
    A[11]=A[10];A[10]=A[9];
    A[9]=A[8];A[8]=A[7];
    A[7]=A[6];A[6]=A[5];
    A[5]=A[4];A[4]=A[3];
    A[3]=A[2];A[2]=A[1];
    A[1]=A[0];
}

void op_rs2(void)
{
    op_rs1();
    op_rs1();
}

void op_srd(void)
{
    set(A,mem[D][to_num(E,6)]);
    op_ls1();
    to_mem(D, E,6, A);
}

void op_sri(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    ad %= 07777;
    
    set(A,mem[I][ad]);
    op_ls1();
    to_mem(I,E,6,A);
}

void op_srf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    set(A,mem[R][ad]);
    op_ls1();
    set(mem[R][ad],A);
}

void op_srb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    set(A,mem[R][ad]);
    op_ls1();
    set(mem[R][ad],A);
}

void op_srs(void)
{
    set(A,mem[0][07777]);
    op_ls1();
    set(mem[0][07777],A);
}

void op_src(void)
{
    set(A,G);
    op_ls1();
    set(G,A);
}

void op_srm(void)
{
    set(A,mem[I][to_num(G,12)]);
    op_ls1();
    to_mem(I,G,12,A);
}

void op_lpn(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for(int i=6;i<RS;i++) {
        t[i] = E[i-6];
    }
    lp(A,t);
}

void op_lpd(void)
{
    lp(A,mem[D][to_num(E,6)]);
}

void op_lpi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    ad %= 07777;
    
    lp(A,mem[I][ad]);
}

void op_lpf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    lp(A,mem[R][ad]);
}

void op_lpb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    lp(A,mem[R][ad]);
}

void op_lps(void)
{
    lp(A,mem[0][07777]);
}

void op_lpc(void)
{
    lp(A,G);
}

void op_lpm(void)
{
    lp(A,mem[I][to_num(G,12)]);
}

void op_scn(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for(int i=6;i<RS;i++) {
        t[i] = E[i-6];
    }
    for(int i=0;i<RS;i++) {
        if (t[i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scd(void)
{
    for(int i=0;i<RS;i++) {
        if(mem[D][to_num(E,6)][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_sci(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    for(int i=0;i<RS;i++) {
        if (mem[I][ad][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scf(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        if(mem[R][ad][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scb(void)
{
    int ad = 0;
    ad = to_num(P,12) - to_num(E,6);
    ad %= 07777;
    
    for(int i=0;i<RS;i++) {
        if (mem[R][ad][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scs(void)
{
    for(int i=0;i<RS;i++) {
        if (mem[0][07777][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scc(void)
{
    for(int i=0;i<RS;i++) {
        if (G[i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scm(void)
{
    for(int i=0;i<RS;i++) {
        if (mem[I][to_num(G,12)][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_srj(void)
{
    R = E[3]*4+E[4]*2+E[5];
    set(P,A);
}

void op_sic(void)
{
    I = E[3]*4+E[4]*2+E[5];
}

void op_irj(void)
{
    I = E[3]*4+E[4]*2+E[5];
    R = E[3]*4+E[4]*2+E[5];
    set(P,A);
}

void op_sdc(void)
{
    D = E[3]*4+E[4]*2+E[5];
}

void op_drj(void)
{
    D = E[3]*4+E[4]*2+E[5];
    R = E[3]*4+E[4]*2+E[5];
    set(P,A);
}

void op_sid(void)
{
    D = E[3]*4+E[4]*2+E[5];
    I = E[3]*4+E[4]*2+E[5];
}

void op_acj(void)
{
    D = E[3]*4+E[4]*2+E[5];
    I = E[3]*4+E[4]*2+E[5];
    R = E[3]*4+E[4]*2+E[5];
    set(P,A);
}

void op_sbu(void)
{
    B = E[3]*4+E[4]*2+E[5];
}

void op_zjf(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(!ne(A,zr0)) {
        adder(P,EE,1);
    } else {
        inc(P);
    }
}

void op_nzf(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(ne(A,zr0)) {
        adder(P,EE,1);
    } else {
        inc(P);
    }
}

void op_pjf(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(A[0] == 0) {
        adder(P,EE,1);
    } else {
        inc(P);
    }
}

void op_njf(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(A[0] == 1) {
        adder(P,EE,1);
    } else {
        inc(P);
    }
}

void op_zjb(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(!ne(A,zr0)) {
        adder(P,EE,0);
    } else {
        inc(P);
    }
}

void op_nzb(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(ne(A,zr0)) {
        adder(P,EE,0);
    } else {
        inc(P);
    }
}

void op_pjb(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(A[0] == 0) {
        adder(P,EE,0);
    } else {
        inc(P);
    }
}

void op_njb(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    cpl(EE,E);
    if(A[0] == 1) {
        adder(P,EE,0);
    } else {
        inc(P);
    }
}

void op_jpi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E,6)],12);
    
    set(P,mem[R][ad]);
}

void op_jpr(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    set(t,P);
    inc(t); inc(t);
    to_mem(R,G,12,t);
    inc(G);
    set(P,G);
}

void op_jfi(void)
{
    int ad = 0;
    ad = to_num(P,12) + to_num(E,6);
    ad %= 07777;
    
    set(P,mem[R][ad]);
}

void op_cbc(void)
{
    // clear buffer controls
    // interrupt any I/O operation
    BFR_busy = 0;
}

void op_cil(void)
{
    // clear interrupt lockout
}

void op_ibi(void)
{
    if (BFR_busy) {
        set(P,G);
    } else {
        int bxr=to_num(BXR,12);
        int ber=to_num(BER,12);
        for(int i=ber;i<=bxr;i++) {
            set(mem[B][i],BFR);
        }
    }
}

void op_ibo(void)
{
    if (BFR_busy) {
        set(P,G);
    } else {
        int bxr=to_num(BXR,12);
        int ber=to_num(BER,12);
        for(int i=ber;i<=bxr;i++) {
            set(BFR,mem[B][i]);
        }
    }
}

void op_inp(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    int fwa = 0;
    int lwa = 0;
    
    cpl(EE,E);
    adder(P,EE,1);
    fwa = to_num(mem[R][to_num(P,12)],12);
    lwa = to_num(G,12);
    for(int i=fwa;i<=lwa;i++) {
        set(mem[I][i],rdd(dvc));
    }
}

void op_out(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    int fwa = 0;
    int lwa = 0;
    
    cpl(EE,E);
    adder(P,EE,1);
    fwa = to_num(mem[R][to_num(P,12)],12);
    lwa = to_num(G,12);
    for(int i=fwa;i<=lwa;i++) {
        wrd(dvc,mem[I][i]);
    }
}

void op_otn(void)
{
    int EE[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    cpl(EE,E);
    wrd(dvc,EE);
}
            
void op_exc(void)
{
    int ret[] = {1,0,0,0,0,0,0,0,0,0,0,0}; // 4000
    
    switch(to_num(G,12)) {
        case 04102:
            dvc = fopen("ptrdr.in","r");
            dvcstatus = 0;
            break;
        case 04104:
            dvc = fopen("ptpch.out","w");
            dvcstatus = 0;
            break;
        case 0600:
            dvc = fopen("prt.out","rw");
            dvcstatus = 1;
            if (dvc == NULL) {
                set(DW,zr0);
            } else {
                set(DW,ret);
            }
            break;
        case 0601:
            fprintf(dvc,"\n");
            break;
        case 0602:
            fprintf(dvc,"\n\n");
            break;
        case 0605:
            fprintf(dvc,"INFO ### A ");
            for(int i=0;i<RS;i++)
                fprintf(dvc,"%d",A[i]);
            fprintf(dvc,"### P ");
            for(int i=0;i<RS;i++)
                fprintf(dvc,"%d",P[i]);
            fprintf(dvc,"\n");
            fprintf(dvc,"\n");
            fprintf(dvc,"\n");
        default:
            break;
    }
}

void op_exf(void)
{
    int ret[] = {1,0,0,0,0,0,0,0,0,0,0,0}; // 4000
    int ad = to_num(P,12);
    
    ad += to_num(E,6);
    
    switch(ad) {
        case 04102:
            dvc = fopen("ptrdr.in","r");
            dvcstatus = 0;
            break;
        case 04104:
            dvc = fopen("ptpch.out","w");
            dvcstatus = 0;
            break;
        case 0600:
        case 0607:
            dvc = fopen("prt.out","rw");
            dvcstatus = 1;
            if (dvc == NULL) {
                set(DW,zr0);
            } else {
                set(DW,ret);
            }
            break;
        case 0601:
            fprintf(dvc,"\n");
            break;
        case 0602:
            fprintf(dvc,"\n\n");
            break;
        case 0605:
            fprintf(dvc,"INFO ### A ");
            for(int i=0;i<RS;i++)
                fprintf(dvc,"%d",A[i]);
            fprintf(dvc,"### P ");
            for(int i=0;i<RS;i++)
                fprintf(dvc,"%d",P[i]);
            fprintf(dvc,"\n");
            fprintf(dvc,"\n");
            fprintf(dvc,"\n");
        default:
            break;
    }
}

void op_sls(void)
{
    int sw = E[3]*4+E[4]*2+E[5];
    
    switch(sw) {
        case 1:
            if (sw1) is_hlt = 1;
            break;
        case 2:
            if (sw2) is_hlt = 1;
            break;
        case 3:
            if (sw1 || sw2) is_hlt = 1;
            break;
        case 4:
            if (sw4) is_hlt = 1;
            break;
        case 5:
            if (sw1 | sw4) is_hlt = 1;
            break;
        case 6:
            if (sw2 | sw4) is_hlt = 1;
            break;
        case 7:
            if (sw1 | sw2 | sw4) is_hlt = 1;
            break;
    }
}

void op_slj(void)
{
    int sw = E[0]*4+E[1]*2+E[2];
    
    switch(sw) {
        case 1:
            if (sw1) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 2:
            if (sw2) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 3:
            if (sw1 || sw2) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 4:
            if (sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 5:
            if (sw1 | sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 6:
            if (sw2 | sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 7:
            if (sw1 | sw2 | sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        default:
            break;
    }
}

void op_sjs(void)
{
    int sw = E[0]*4+E[1]*2+E[2];
    
    switch(sw) {
        case 1:
            if (sw1) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 2:
            if (sw2) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 3:
            if (sw1 || sw2) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 4:
            if (sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 5:
            if (sw1 | sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 6:
            if (sw2 | sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        case 7:
            if (sw1 | sw2 | sw4) {
                set(P,G);
            } else {
                inc(P);
            }
            break;
        default:
            break;
    }
    
    sw = E[3]*4+E[4]*2+E[5];
    
    switch(sw) {
        case 1:
            if (sw1) is_hlt = 1;
            break;
        case 2:
            if (sw2) is_hlt = 1;
            break;
        case 3:
            if (sw1 || sw2) is_hlt = 1;
            break;
        case 4:
            if (sw4) is_hlt = 1;
            break;
        case 5:
            if (sw1 | sw4) is_hlt = 1;
            break;
        case 6:
            if (sw2 | sw4) is_hlt = 1;
            break;
        case 7:
            if (sw1 | sw2 | sw4) is_hlt = 1;
            break;
    }
}

void microcode(void)
{
    readFE();
    switch(to_num(F,6)) {
        case 0:
            switch(to_num(E,6)) {
                case 0:
                    op_err();
                    break;
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                    op_nop();
                    inc(P);
                    break;
                case 010:
                case 011:
                case 012:
                case 013:
                case 014:
                case 015:
                case 016:
                case 017:
                    op_srj();
                    break;
                case 020:
                case 021:
                case 022:
                case 023:
                case 024:
                case 025:
                case 026:
                case 027:
                    op_sic();
                    inc(P);
                    break;
                case 030:
                case 031:
                case 032:
                case 033:
                case 034:
                case 035:
                case 036:
                case 037:
                    op_irj();
                    break;
                case 040:
                case 041:
                case 042:
                case 043:
                case 044:
                case 045:
                case 046:
                case 047:
                    op_sdc();
                    inc(P);
                    break;
                case 050:
                case 051:
                case 052:
                case 053:
                case 054:
                case 055:
                case 056:
                case 057:
                    op_drj();
                    break;
                case 060:
                case 061:
                case 062:
                case 063:
                case 064:
                case 065:
                case 066:
                case 067:
                    op_sid();
                    inc(P);
                    break;
                case 070:
                case 071:
                case 072:
                case 073:
                case 074:
                case 075:
                case 076:
                case 077:
                    op_acj();
                    break;
            }
            break;
        case 1:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_bls();
                    inc(P);
                    break;
                case 1:
                    op_pta();
                    inc(P);
                    break;
                case 2:
                    op_ls1();
                    inc(P);
                    break;
                case 3:
                    op_ls2();
                    inc(P);
                    break;
                case 5:
                    readG();
                    op_ate();
                    inc(P);
                    break;
                case 6:
                    readG();
                    op_atx();
                    inc(P);
                    break;
                case 7:
                    op_eta();
                    inc(P);
                    break;
                case 010:
                    op_ls3();
                    inc(P);
                    break;
                case 011:
                    op_ls6();
                    inc(P);
                    break;
                case 012:
                    op_mut();
                    inc(P);
                    break;
                case 013:
                    op_muh();
                    inc(P);
                    break;
                case 014:
                    op_rs1();
                    inc(P);
                    break;
                case 015:
                    op_rs2();
                    inc(P);
                    break;
                case 030:
                    op_cta();
                    inc(P);
                    break;
                case 040:
                case 041:
                case 042:
                case 043:
                case 044:
                case 045:
                case 046:
                case 047:
                    op_sbu();
                    inc(P);
                    break;
                case 050:
                case 051:
                case 052:
                case 053:
                case 054:
                case 055:
                case 056:
                case 057:
                    op_stp();
                    inc(P);
                    break;
                case 060:
                case 061:
                case 062:
                case 063:
                case 064:
                case 065:
                case 066:
                case 067:
                    op_ste();
                    inc(P);
                    break;
            }
            break;
        case 2:
            op_lpn();
            inc(P);
            break;
        case 3:
            op_scn();
            inc(P);
            break;
        case 4:
            op_ldn();
            inc(P);
            break;
        case 5:
            op_lcn();
            inc(P);
            break;
        case 6:
            op_adn();
            inc(P);
            break;
        case 7:
            op_sbn();
            inc(P);
            break;
        case 010:
            op_lpd();
            inc(P);
            break;
        case 011:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_lpm();
                    inc(P);
                    break;
                default:
                    op_lpi();
                    inc(P);
                    break;
            }
            break;
        case 012:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_lpc();
                    inc(P);
                    break;
                default:
                    op_lpf();
                    inc(P);
                    break;
            }
            break;
        case 013:
            switch(to_num(E,6)) {
                case 0:
                    op_lps();
                    inc(P);
                    break;
                default:
                    op_lpb();
                    inc(P);
                    break;
            }
            break;
        case 014:
            op_scd();
            inc(P);
            break;
        case 015:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_scm();
                    inc(P);
                    break;
                default:
                    op_sci();
                    inc(P);
                    break;
            }
            break;
        case 016:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_scc();
                    inc(P);
                    break;
                default:
                    op_scf();
                    inc(P);
                    break;
            }
            break;
        case 017:
            switch(to_num(E,6)) {
                case 0:
                    op_scs();
                    inc(P);
                    break;
                default:
                    op_scb();
                    inc(P);
                    break;
            }
            break;
        case 020:
            op_ldd();
            inc(P);
            break;
        case 021:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_ldm();
                    inc(P);
                    break;
                default:
                    op_ldi();
                    inc(P);
                    break;
            }
            break;
        case 022:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_ldc();
                    inc(P);
                    break;
                default:
                    op_ldf();
                    inc(P);
                    break;
            }
            break;
        case 023:
            switch(to_num(E,6)) {
                case 0:
                    op_lds();
                    inc(P);
                    break;
                default:
                    op_ldb();
                    inc(P);
                    break;
            }
            break;
        case 024:
            op_lcd();
            inc(P);
            break;
        case 025:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_lcm();
                    inc(P);
                    break;
                default:
                    op_lci();
                    inc(P);
                    break;
            }
            break;
        case 026:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_lcc();
                    inc(P);
                    break;
                default:
                    op_lcf();
                    inc(P);
                    break;
            }
            break;
        case 027:
            switch(to_num(E,6)) {
                case 0:
                    op_lcs();
                    inc(P);
                    break;
                default:
                    op_lcb();
                    inc(P);
                    break;
            }
            break;
        case 030:
            op_add();
            inc(P);
            break;
        case 031:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_adm();
                    inc(P);
                    break;
                default:
                    op_adi();
                    inc(P);
                    break;
            }
            break;
        case 032:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_adc();
                    inc(P);
                    break;
                default:
                    op_adf();
                    inc(P);
                    break;
            }
            break;
        case 033:
            switch(to_num(E,6)) {
                case 0:
                    op_ads();
                    inc(P);
                    break;
                default:
                    op_adb();
                    inc(P);
                    break;
            }
            break;
        case 034:
            op_sbd();
            inc(P);
            break;
        case 035:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_sbm();
                    inc(P);
                    break;
                default:
                    op_sbi();
                    inc(P);
                    break;
            }
            break;
        case 036:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_sbc();
                    inc(P);
                    break;
                default:
                    op_sbf();
                    inc(P);
                    break;
            }
            break;
        case 037:
            switch(to_num(E,6)) {
                case 0:
                    op_sbs();
                    inc(P);
                    break;
                default:
                    op_sbb();
                    inc(P);
                    break;
            }
            break;
        case 040:
            op_std();
            inc(P);
            break;
        case 041:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_stm();
                    inc(P);
                    break;
                default:
                    op_sti();
                    inc(P);
                    break;
            }
            break;
        case 042:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_stc();
                    inc(P);
                    break;
                default:
                    op_stf();
                    inc(P);
                    break;
            }
            break;
        case 043:
            switch(to_num(E,6)) {
                case 0:
                    op_sts();
                    inc(P);
                    break;
                default:
                    op_stb();
                    inc(P);
                    break;
            }
            break;
        case 044:
            op_srd();
            inc(P);
            break;
        case 045:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_srm();
                    inc(P);
                    break;
                default:
                    op_sri();
                    inc(P);
                    break;
            }
            break;
        case 046:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_src();
                    inc(P);
                    break;
                default:
                    op_srf();
                    inc(P);
                    break;
            }
            break;
        case 047:
            switch(to_num(E,6)) {
                case 0:
                    op_srs();
                    inc(P);
                    break;
                default:
                    op_srb();
                    inc(P);
                    break;
            }
            break;
        case 050:
            op_rad();
            inc(P);
            break;
        case 051:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_ram();
                    inc(P);
                    break;
                default:
                    op_rai();
                    inc(P);
                    break;
            }
            break;
        case 052:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_rac();
                    inc(P);
                    break;
                default:
                    op_raf();
                    inc(P);
                    break;
            }
            break;
        case 053:
            switch(to_num(E,6)) {
                case 0:
                    op_ras();
                    inc(P);
                    break;
                default:
                    op_rab();
                    inc(P);
                    break;
            }
            break;
        case 054:
            op_aod();
            inc(P);
            break;
        case 055:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_aom();
                    inc(P);
                    break;
                default:
                    op_aoi();
                    inc(P);
                    break;
            }
            break;
        case 056:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_aoc();
                    inc(P);
                    break;
                default:
                    op_aof();
                    inc(P);
                    break;
            }
            break;
        case 057:
            switch(to_num(E,6)) {
                case 0:
                    op_aos();
                    inc(P);
                    break;
                default:
                    op_aob();
                    inc(P);
                    break;
            }
            break;
        case 060:
            op_zjf();
            break;
        case 061:
            op_nzf();
            break;
        case 062:
            op_pjf();
            break;
        case 063:
            op_njf();
            break;
        case 064:
            op_zjb();
            break;
        case 065:
            op_nzb();
            break;
        case 066:
            op_pjb();
            break;
        case 067:
            op_njb();
            break;
        case 070:
            op_jpi();
            break;
        case 071:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_jpr();
                    break;
                default:
                    op_jfi();
                    break;
            }
            break;
        case 072:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_ibi();
                    inc(P);
                    break;
                default:
                    readG();
                    op_inp();
                    inc(P);
                    break;
            }
            break;
        case 073:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_ibo();
                    inc(P);
                    break;
                default:
                    readG();
                    op_out();
                    inc(P);
                    break;
            }
            break;
        case 074:
            op_otn();
            inc(P);
            break;
        case 075:
            switch(to_num(E,6)) {
                case 0:
                    readG();
                    op_exc();
                    inc(P);
                    break;
                default:
                    op_exf();
                    inc(P);
                    break;
            }
            break;
        case 076:
            switch(to_num(E,6)) {
                case 0:
                    op_ina();
                    inc(P);
                    break;
                case 077:
                    op_ota();
                    inc(P);
                    break;
                default:
                    op_hwi();
                    inc(P);
            }
            break;
        case 077:
            switch(to_num(E,6)) {
                case 0:
                    op_hlt(0);
                    break;
                case 077:
                    op_hlt(077);
                    break;
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                case 6:
                case 7:
                    op_sls();
                    inc(P);
                    break;
                case 010:
                case 020:
                case 030:
                case 040:
                case 050:
                case 060:
                case 070:
                    readG();
                    op_slj();
                    break;
                default:
                    readG();
                    op_sjs();
                    break;
            }
            break;
    }
}

void _dmem(int bl, int s, int e)
{
    printf("Memory bank %d dump\n",bl);
    for(int i=s,pos=1;i<=e;i++,pos++) {
        printf("[%04o]  ",i);
        for(int j=0;j<12;j++) {
            printf("%d",mem[bl][i][j]);
        }
        if (pos%4 == 0) {
            printf("\n");
        } else {
            printf(" ");
        }
    }
    printf(".\n");
}

void dump(void)
{
    printf("\n\nDUMPING MEMORY INFO\n");
    printf("------- ------ ----\n");
    printf("A  ");
    for(int i=0;i<RS;i++)
        printf("%d",A[i]);
    printf(" (%04o)",to_num(A,12));
    printf("\n");
    printf("P ");
    printf("%d",R);
    for(int i=0;i<RS;i++)
        printf("%d",P[i]);
    printf(" (%04o)",to_num(P,12));
    printf("\n");
    printf("S  ");
    for(int i=0;i<RS;i++)
        printf("%d",mem[0][07777][i]);
    printf("\n");
    _dmem(0,0,0700);
    //_dmem(1,0100,07771);
}

void ldmp(const char* fn)
{
    FILE* f=fopen(fn,"r");
    char c;
    int i = 0;
    int j = 0;
    int org = 0;
    int bl = 0;
    int BOP = 0;
    
    if (f == NULL) {
        printf("Cannot open memory file!\n");
        exit(100);
    }
    
    printf("Opening %s...\n",fn);
    printf("Memory locations used (octal)\n");
    
    while((c=fgetc(f)) != EOF) {
        if (c == '#') {
            bl = fgetc(f)-'0';
            fgetc(f);
            org = (fgetc(f)-'0') * 512 + (fgetc(f)-'0') * 64 +
                    (fgetc(f)-'0') * 8 + (fgetc(f)-'0');
            printf("%04o\n",org);
            j = 0;
        }
        else if (c == '@') {
            BOP = (fgetc(f)-'0') * 512 + (fgetc(f)-'0') * 64 +
            (fgetc(f)-'0') * 8 + (fgetc(f)-'0');
        }
        else if ((c != ' ') && (c != '\n')) {
            mem[bl][org+j][i] = c - '0';
            i++;
            if ((i % 12) == 0) {
                i = 0;
                j++;
                printf("%04o ",org+j);
                if (j%6 == 0) {
                    printf("\n");
                }
            }
        }
    }
    _set(P,BOP);
    printf("\nP set to %04o",to_num(P,12));
    printf("\nLoading completed.\n\n");
}

int main(int argc, const char * argv[]) {
    init();
    
    if (argc > 1) {
        printf("Loading memory dmp...\n");
        ldmp(argv[1]);
    }
    
    
    //_set(A,0100);
    
    while (is_hlt == 0) {
        printf("A  ");
        for(int i=0;i<RS;i++)
            printf("%d",A[i]);
        printf(" (%04o)",to_num(A,12));
        printf("\n");
        printf("P ");
        printf("%d",R);
        for(int i=0;i<RS;i++)
            printf("%d",P[i]);
        printf(" (%04o)",to_num(P,12));
        printf("\n---\n");
        microcode();
    }
    
    fclose(tape);
    fclose(prt);
    
    return 0;
}

