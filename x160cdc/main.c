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

int* sbc;

int*** mem;

#define RS 12

int one[] = {0,0,0,0,0,0,0,0,0,0,0,1};

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

void adder(int* a, int* b, int op){
    int bi = 0;
    
    if (op)
        cp(b);
    for(int i=RS-1;i>=0;i--) {
        a[i] = fs(a[i],b[i],bi,&ovr);
        bi = ovr;
    }
    if (ovr == 1) {
        bi = 0;
        for(int i=RS-1;i>=0;i--) {
            a[i] = fs(a[i],one[i],bi,&ovr);
        }
    }
}

void inc(int* a){
    int bi = 0;

    cp(one);
    for(int i=RS-1;i>=0;i--) {
        a[i] = fs(a[i],one[i],bi,&ovr);
        bi = ovr;
    }
    if (ovr == 1) {
        bi = 0;
        for(int i=RS-1;i>=0;i--) {
            a[i] = fs(a[i],one[i],bi,&ovr);
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
    A = calloc(12,sizeof(int));
    P = calloc(12,sizeof(int));
    BER = calloc(12,sizeof(int));
    BXR = calloc(12,sizeof(int));
    
    BFR = calloc(12,sizeof(int));
    S = calloc(12,sizeof(int));
    F = calloc(6,sizeof(int));
    E = calloc(6,sizeof(int));
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
}

int to_num(int* b)
{
    int r = 0;
    
    for(int i=11,j=0;i>=0;i--,j++) {
        r += b[i] * pow(2,j);
    }
    
    return r;
}

void set(int* d, int* s)
{
    for(int i=0;i<RS;i++)
        d[i] = s[i];
}

void from_mem(int bl, int* d, int* s)
{
    for(int i=0;i<RS;i++) {
        d[i] = mem[bl][to_num(s)][i];
    }
}

void to_mem(int bl, int* d, int* s)
{
    for(int i=0;i<RS;i++) {
        mem[bl][to_num(d)][i]=s[i];
    }
}

void readG(void)
{
    inc(P);
    from_mem(R,G,P);
}

void op_err(void)
{
    is_hlt = 1;
    printf("ERR");
}

void op_hlt(void)
{
    is_hlt = 1;
}

void op_nop(void)
{
    // nothing, nic!
    inc(P);
}

void op_bls(void)
{
    if (BFR_busy) {
        set(P,G);
    } else {
        BFR_busy = 1;
        set(BFR,A);
        while (ne(BER,BXR)) {
            to_mem(B,BER,BFR);
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
    inc(P);
}

void op_cta(void)
{
    int _B=B;int _D=D;int _I=I;int _R=R;
    
    A[11]=_B%2;_B/=2;A[10]=_B%2;_B/=2;A[9]=_B%2;
    A[8]=_D%2;_D/=2;A[7]=_D%2;_D/=2;A[6]=_D%2;
    A[5]=_I%2;_I/=2;A[4]=_I%2;_I/=2;A[3]=_I%2;
    A[2]=_R%2;_R/=2;A[1]=_R%2;_R/=2;A[0]=_R%2;
    inc(P);
}

void op_stp(void)
{
    for(int i=0;i<11;i++) {
        mem[D][to_num(E)][i] = P[i];
    }
}

void op_ldn(void)
{
    for(int i=0;i<6;i++) {
        A[i] = 0;
    }
    for(int i=6;i<11;i++) {
        A[i] = E[i];
    }
}

void op_ldd(void)
{
    for(int i=0;i<11;i++) {
        A[i] = mem[D][to_num(E)][i];
    }
}

void op_ldi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    for(int i=0;i<11;i++) {
        A[i] = mem[I][ad][i];
    }
}

void op_ldf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    for(int i=0;i<11;i++) {
        A[i] = mem[R][ad][i];
    }
}

void op_ldb(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
    ad %= 07777;
    
    for(int i=0;i<11;i++) {
        A[i] = mem[R][ad][i];
    }
}

void op_lds(void)
{
    for(int i=0;i<11;i++) {
        A[i] = mem[0][07777][i];
    }
}

void op_ldc(void)
{
    set(A,G);
}

void op_ldm(void)
{
    set(A,mem[I][to_num(G)]);
}

void op_lcn(void)
{
    for(int i=0;i<6;i++) {
        A[i] = 0;
    }
    for(int i=6;i<11;i++) {
        A[i] = !E[i];
    }
}

void op_lcd(void)
{
    for(int i=0;i<11;i++) {
        A[i] = !mem[D][to_num(E)][i];
    }
}

void op_lci(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    for(int i=0;i<11;i++) {
        A[i] = !mem[I][ad][i];
    }
}

void op_lcf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    for(int i=0;i<11;i++) {
        A[i] = !mem[R][ad][i];
    }
}
    
void op_lcb(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
    ad %= 07777;
    
    for(int i=0;i<11;i++) {
        A[i] = !mem[R][ad][i];
    }
}

void op_lcs(void)
{
    for(int i=0;i<11;i++) {
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
    set(A,mem[I][to_num(G)]);
    cp(A);
}

void op_ste(void)
{
    to_mem(D,E,BER);
    set(BER,A);
}

void op_std(void)
{
    for(int i=0;i<11;i++) {
        mem[D][to_num(E)][i] = A[i];
    }
}

void op_sti(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    for(int i=0;i<11;i++) {
        mem[I][ad][i] = A[i];
    }
}

void op_stf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    for(int i=0;i<11;i++) {
        mem[R][ad][i] = A[i];
    }
}

void op_stb(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
    ad %= 07777;
    
    for(int i=0;i<11;i++) {
        mem[R][ad][i] = A[i];
    }
}

void op_sts(void)
{
    for(int i=0;i<11;i++) {
        mem[0][07777][i] = A[i];
    }
}

void op_stc(void)
{
    set(G,A);
}

void op_stm(void)
{
    set(mem[I][to_num(G)],A);
}

void op_hwi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    for(int i=6;i<11;i++) {
        mem[I][ad][i] = A[i];
    }
}

void op_ina(void)
{
    
}

void op_ota(void)
{
    
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
    
    for(int i=6;i<11;i++) {
        t[i] = E[i];
    }
    adder(A,t,1);
}

void op_add(void)
{
    adder(A,mem[D][to_num(E)],1);
}

void op_adi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    adder(A,mem[I][ad],1);
}

void op_adf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    adder(A,mem[R][ad],1);
}

void op_adb(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
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
    adder(A,mem[I][to_num(G)],1);
}

void op_sbn(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for(int i=6;i<11;i++) {
        t[i] = E[i];
    }
    adder(A,t,0);
}

void op_sbd(void)
{
    adder(A,mem[D][to_num(E)],0);
}

void op_sbi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    adder(A,mem[I][ad],0);
}

void op_sbf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    adder(A,mem[R][ad],0);
}

void op_sbb(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
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
    adder(A,mem[I][to_num(G)],0);
}

void op_rad(void)
{
    adder(A,mem[D][to_num(E)],1);
    to_mem(D, E, A);
}

void op_rai(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    adder(A,mem[I][ad],1);
    to_mem(I,E,A);
}

void op_raf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    adder(A,mem[R][ad],1);
    set(mem[R][ad],A);
}

void op_rab(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
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
    adder(A,mem[I][to_num(G)],1);
    to_mem(I,G,A);
}

void op_aod(void)
{
    set(A,mem[D][to_num(E)]);
    inc(A);
    to_mem(D, E, A);
}

void op_aoi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    set(A,mem[I][ad]);
    inc(A);
    to_mem(I,E,A);
}

void op_aof(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    set(A,mem[R][ad]);
    inc(A);
    set(mem[R][ad],A);
}

void op_aob(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
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
    set(A,mem[I][to_num(G)]);
    inc(A);
    to_mem(I,G,A);
}

void op_ls1(void)
{
    int t = A[0];
    
    A[0]= A[1]; A[1]= A[2];
    A[3]= A[4]; A[4]= A[5];
    A[6]= A[7]; A[7]= A[8];
    A[9]=A[10];A[10]=A[11];
    A[11]=t;
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
    set(A,mem[D][to_num(E)]);
    op_ls1();
    to_mem(D, E, A);
}

void op_sri(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    set(A,mem[I][ad]);
    op_ls1();
    to_mem(I,E,A);
}

void op_srf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    set(A,mem[R][ad]);
    op_ls1();
    set(mem[R][ad],A);
}

void op_srb(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
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
    set(A,mem[I][to_num(G)]);
    op_ls1();
    to_mem(I,G,A);
}

void op_lpn(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for(int i=6;i<11;i++) {
        t[i] = E[i];
    }
    lp(A,t);
}

void op_lpd(void)
{
    lp(A,mem[D][to_num(E)]);
}

void op_lpi(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    lp(A,mem[I][ad]);
}

void op_lpf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
    ad %= 07777;
    
    lp(A,mem[R][ad]);
}

void op_lpb(void)
{
    int ad = 0;
    ad = to_num(P) - to_num(E);
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
    lp(A,mem[I][to_num(G)]);
}

void op_scn(void)
{
    int t[] = {0,0,0,0,0,0,0,0,0,0,0,0};
    
    for(int i=6;i<RS;i++) {
        t[i] = E[i];
    }
    for(int i=0;i<RS;i++) {
        if (E[i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scd(void)
{
    for(int i=0;i<RS;i++) {
        if(mem[D][to_num(E)][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_sci(void)
{
    int ad = 0;
    ad = to_num(mem[D][to_num(E)]);
    
    for(int i=0;i<RS;i++) {
        if (mem[I][ad][i] == 1) {
            A[i] = !A[i];
        }
    }
}

void op_scf(void)
{
    int ad = 0;
    ad = to_num(P) + to_num(E);
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
    ad = to_num(P) - to_num(E);
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
        if (mem[I][to_num(G)][i] == 1) {
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

void microcode(void)
{
    switch(to_num(F)) {
        case 0:
            switch(to_num(E)) {
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
            }
            break;
        case 1:
            switch(to_num(E)) {
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
            switch(to_num(E)) {
                case 0:
                    readG();
                    op_lpm();
                    inc(P);
                    break;
                default:
                    op_lpi();
                    inc(P);
            }
            break;
        case 012:
            switch(to_num(E)) {
                case 0:
                    readG();
                    op_lpc();
                    inc(P);
                    break;
                default:
                    op_lpf();
                    inc(P);
            }
            break;
        case 013:
            switch(to_num(E)) {
                case 0:
                    op_lps();
                    inc(P);
                    break;
                default:
                    op_lpb();
                    inc(P);
            }
            break;
        case 014:
            op_scd();
            inc(P);
            break;
        case 015:
            switch(to_num(E)) {
                case 0:
                    readG();
                    op_scm();
                    inc(P);
                    break;
                default:
                    op_sci();
                    inc(P);
            }
            break;
        case 016:
            switch(to_num(E)) {
                case 0:
                    readG();
                    op_scc();
                    inc(P);
                    break;
                default:
                    op_scf();
                    inc(P);
            }
            break;
        case 017:
            switch(to_num(E)) {
                case 0:
                    op_scs();
                    inc(P);
                    break;
                default:
                    op_scb();
                    inc(P);
            }
            break;
        case 020:
            op_ldd();
            inc(P);
            break;
        case 021:
            switch(to_num(E)) {
                case 0:
                    readG();
                    op_ldm();
                    inc(P);
                    break;
                default:
                    op_ldi();
                    inc(P);
            }
            break;
        case 022:
            switch(to_num(E)) {
                case 0:
                    readG();
                    op_ldc();
                    inc(P);
                    break;
                default:
                    op_ldf();
                    inc(P);
            }
            break;
        case 023:
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
            switch(to_num(E)) {
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
        case 076:
            switch(to_num(E)) {
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
        case 077:
            switch(to_num(E)) {
                case 0:
                case 077:
                    op_hlt();
                    break;
            }
            break;
    }
}

void dump(void)
{
    for(int i=0;i<RS;i++)
        printf("%d",A[i]);
}

int main(int argc, const char * argv[]) {
    // insert code here...
    init();
    
    // 000 000 110 011 (0063)
    A[0] = 0;  A[1] = 0;  A[2] = 0;
    A[3] = 0;  A[4] = 0;  A[5] = 0;
    A[6] = 1;  A[7] = 1;  A[8] = 0;
    A[9] = 0; A[10] = 1; A[11] = 1;
    
    // 000 000 010 011 (0023)
    Z[0] = 0;  Z[1] = 0;  Z[2] = 0;
    Z[3] = 0;  Z[4] = 0;  Z[5] = 0;
    Z[6] = 0;  Z[7] = 1;  Z[8] = 0;
    Z[9] = 0; Z[10] = 1; Z[11] = 1;
    
    adder(A,Z,0);
    
    dump();
    
    return 0;
}
