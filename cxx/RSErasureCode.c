/*
Encoding/erasure decoding for Reed-Solomon codes over binary extension fields
Author: Sian-Jheng Lin (King Abdullah University of Science and Technology (KAUST), email: sianjheng.lin@kaust.edu.sa)

This program is the implementation of
Lin, Han and Chung, "Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes," FOCS14.
(http://arxiv.org/abs/1404.3458)
*/

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdio.h>


#include "sha-256.h"

/*
typedef unsigned char GFSymbol;
#define len 8//2^len: the size of Galois field
GFSymbol mask = 0x1D; //GF(2^8): x^8 + x^4 + x^3 + x^2 + 1
GFSymbol Base[] = {1, 214, 152, 146, 86, 200, 88, 230};//Cantor basis
*/

typedef unsigned short GFSymbol;
#define len 16
GFSymbol mask = 0x2D;//x^16 + x^5 + x^3 + x^2 + 1
GFSymbol Base[len] = {1, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198};//Cantor basis

#define Size (1<<len)//Field size
#define mod (Size-1)

GFSymbol log_tbl[Size];
GFSymbol exp_tbl[Size];

//-----Used in decoding procedure-------
GFSymbol skewVec[mod];//twisted factors used in FFT
GFSymbol B[Size>>1];//factors used in formal derivative
GFSymbol log_walsh[Size];//factors used in the evaluation of the error locator polynomial

GFSymbol mulE(GFSymbol a, GFSymbol b){//return a*exp_tbl[b] over GF(2^r)
	return a? exp_tbl[(log_tbl[a]+b &mod) + (log_tbl[a]+b >>len)]: 0;
}

void walsh(GFSymbol* data, int size){//fast Walsh–Hadamard transform over modulo mod
	for (int depart_no=1; depart_no<size; depart_no <<= 1){
		for (int j = 0; j < size; j += depart_no<<1){
			for (int i=j; i<depart_no+j; i++){
				unsigned tmp2 = data[i] + mod - data[i+depart_no];
				data[i] = (data[i] + data[i+depart_no]&mod) + (data[i] + data[i+depart_no]>>len);
				data[i+depart_no] = (tmp2&mod) + (tmp2>>len);
			}
		}
	}
	return;
}

void formal_derivative(GFSymbol* cos, int size){//formal derivative of polynomial in the new basis
	for(int i=1; i<size; i++){
		int leng = ((i^i-1)+1)>>1;
		for(int j=i-leng; j<i; j++)
			cos[j] ^= cos[j+leng];
	}
	for(int i=size; i<Size; i<<=1)
		for(int j=0; j<size; j++)
			cos[j] ^= cos[j+i];
	return;
}

void IFLT(GFSymbol* data, int size, int index){//IFFT in the proposed basis
	for (int depart_no=1; depart_no<size; depart_no <<= 1){
		for (int j=depart_no; j < size; j += depart_no<<1){
			for (int i=j-depart_no; i<j; i++)
				data[i+depart_no] ^= data[i];

			GFSymbol skew = skewVec[j+index-1];
			if (skew != mod)
				for (int i=j-depart_no; i<j; i++)
					data[i] ^= mulE(data[i+depart_no], skew);
		}
	}
	return;
}

void FLT(GFSymbol* data, int size, int index){//FFT in the proposed basis
	for(int depart_no = size>>1; depart_no > 0; depart_no >>= 1){
		for (int j = depart_no; j < size; j += depart_no<<1){
			GFSymbol skew = skewVec[j+index-1];
			if (skew != mod)
				for (int i=j-depart_no; i<j; i++)
					data[i] ^= mulE(data[i+depart_no], skew);
			for (int i=j-depart_no; i<j; i++)
				data[i+depart_no] ^= data[i];
		}
	}
	return;
}

void init(){//initialize log_tbl[], exp_tbl[]
	GFSymbol mas = (1<<len-1)-1;
	GFSymbol state=1;
	for(int i=0; i<mod; i++){
		exp_tbl[state]=i;
        if(state>>len-1){
        	state &= mas;
        	state = state<<1^mask;
        }else
        	state <<= 1;
    }
    exp_tbl[0] = mod;

    log_tbl[0] = 0;
	for(int i=0; i<len; i++)
		for(int j=0; j<1<<i; j++)
			log_tbl[j+(1<<i)] = log_tbl[j] ^ Base[i];
    for(int i=0; i<Size; i++)
        log_tbl[i]=exp_tbl[log_tbl[i]];

    for(int i=0; i<Size; i++)
        exp_tbl[log_tbl[i]]=i;
    exp_tbl[mod] = exp_tbl[0];
}


void init_dec(){//initialize skewVec[], B[], log_walsh[]
	GFSymbol base[len-1];

	for(int i=1; i<len; i++)
		base[i-1] = 1<<i;

	for(int m=0; m<len-1; m++){
		int step = 1<<(m+1);
		skewVec[(1<<m)-1] = 0;
		for(int i=m; i<len-1; i++){
			int s = 1<<(i+1);
			for(int j=(1<<m)-1; j<s; j+=step)
				skewVec[j+s] = skewVec[j] ^ base[i];
		}
		base[m] = mod-log_tbl[mulE(base[m], log_tbl[base[m]^1])];
		for(int i=m+1; i<len-1; i++)
			base[i] = mulE(base[i], (log_tbl[base[i]^1]+base[m])%mod);
	}
	for(int i=0; i<Size; i++)
		skewVec[i] = log_tbl[skewVec[i]];

	base[0] = mod-base[0];
	for(int i=1; i<len-1; i++)
		base[i] = (mod-base[i]+base[i-1])%mod;

	B[0] = 0;
	for(int i=0; i<len-1; i++){
		int depart = 1<<i;
		for(int j=0; j<depart; j++)
			B[j+depart] = (B[j] + base[i])%mod;
	}

	memcpy(log_walsh, log_tbl, Size*sizeof(GFSymbol));
	log_walsh[0] = 0;
	walsh(log_walsh, Size);
}

void encodeL(GFSymbol* data, int k, GFSymbol* codeword){//Encoding alg for k/n<0.5: message is a power of two
	memcpy(codeword, data, sizeof(GFSymbol)*k);
	IFLT(codeword, k, 0);
	for(int i=k; i<Size; i+=k){
		memcpy(&codeword[i], codeword, sizeof(GFSymbol)*k);
		FLT(&codeword[i], k, i);
	}
	memcpy(codeword, data, sizeof(GFSymbol)*k);
	return;
}

void encodeH(GFSymbol* data, int k, GFSymbol* parity, GFSymbol* mem){//Encoding alg for k/n>0.5: parity is a power of two.
//data: message array. parity: parity array. mem: buffer(size>= n-k)
	int t = Size-k;
	memset(parity, 0, sizeof(GFSymbol)*t);
	for(int i=t; i<Size; i+=t){
		memcpy(mem, &data[i-t], sizeof(GFSymbol)*t);
		IFLT(mem, t, i);
		for(int j=0; j<t; j++)
			parity[j] ^= mem[j];
	}
	FLT(parity, t, 0);
	return;
}

//Compute the evaluations of the error locator polynomial
void decode_init(_Bool* erasure, GFSymbol* log_walsh2){
	for(int i=0; i<Size; i++)
		log_walsh2[i] = erasure[i];
	walsh(log_walsh2, Size);
	for (int i=0; i<Size; i++)
		log_walsh2[i] = (unsigned long)log_walsh2[i]*log_walsh[i]%mod;
	walsh(log_walsh2,Size);
	for (int i=0; i<Size; i++)
		if(erasure[i]) log_walsh2[i] = mod-log_walsh2[i];
}

void decode_main(GFSymbol* codeword, _Bool* erasure, GFSymbol* log_walsh2){
	int k2 = Size;//k2 can be replaced with k
	for (int i=0; i<Size; i++)
		codeword[i] = erasure[i]? 0 : mulE(codeword[i], log_walsh2[i]);
	IFLT(codeword, Size, 0);

	for(int i=0; i<Size; i+=2) {//formal derivative
		codeword[i] = mulE(codeword[i], mod-B[i>>1]);
		codeword[i+1] = mulE(codeword[i+1], mod-B[i>>1]);
	}

	formal_derivative(codeword, k2);
	for(int i=0; i<k2; i+=2){
		codeword[i] = mulE(codeword[i], B[i>>1]);
		codeword[i+1] = mulE(codeword[i+1], B[i>>1]);
	}

	FLT(codeword, k2, 0);
	for (int i=0; i<k2; i++) {
		codeword[i] = erasure[i]? mulE(codeword[i], log_walsh2[i]) : 0;
	}
}

const int N = 1<<8;


void print_sha256(char* txt, uint8_t* data, size_t lx) {
	uint8_t hash[32] = {0, };
	calc_sha_256(hash, data, lx);
	printf("sha256(c|%s):\n", txt);
	for(int i=0; i<32; i++) {
		printf("%02x", hash[i]);
	}
	printf("\n");
}

void test(int k){
	//-----------Generating message----------
	GFSymbol data[Size] = {0};//message array
	srand(time(NULL));
	// for(int i=Size-k; i<Size; i++)
    for(int i=0; i<k; i++)
		data[i] = i*i % mod; // rand()&mod;//filled with random numbers


	printf("Message(Last n-k are zeros): \n");
	for(int i=0; i<k; i++) {
		printf("%04x ", data[i]);
	}
	printf("\n");

	print_sha256("data", (uint8_t*)data, Size * 2);

	//---------encoding----------
	GFSymbol codeword[Size];
	// encodeH(&data[Size-k], k, &data, codeword);
	encodeL(data, k, codeword);

	print_sha256("encoded", (uint8_t*)codeword, Size*2);

	// memcpy(codeword, data, sizeof(GFSymbol)*Size);

	printf("Codeword:\n");
	for(int i=k; i<(k+100); i++) {
		printf("%02X ", codeword[i]);
	}
	printf("\n");

	//--------erasure simulation---------
	_Bool erasure[Size] = {0};//Array indicating erasures
	for(int i=0; i<(Size-k); i++)
		erasure[i] = 1;

	#if 0
		for(int i=Size-1; i>0; i--){//permuting the erasure array
			int pos = rand()%(i+1);
			if(i != pos){
				_Bool tmp = erasure[i];
				erasure[i] = erasure[pos];
				erasure[pos] = tmp;
			}
		}
	#endif

	for (int i=0; i<Size; i++) {
		if(erasure[i]) codeword[i] = 0;
	}


	print_sha256("erased", (uint8_t*)codeword, Size*2);

	printf("Erasure (XX is erasure):\n");
	for(int i=0; i<Size; i++) {
		// if(erasure[i]) printf("XX ");
		// else printf("%02X ", codeword[i]);
	}
	printf("\n");


	//---------Erasure decoding----------------
	GFSymbol log_walsh2[Size];
	decode_init(erasure, log_walsh2);//Evaluate error locator polynomial
	print_sha256("log_walsh2", (uint8_t*)log_walsh2, Size*2);

	//---------main processing----------
	decode_main(codeword, erasure, log_walsh2);

	print_sha256("recovered", (uint8_t*)codeword, Size*2);

	printf("Decoded result:\n");
	for(int i=0; i<(k+10); i++){
		if(erasure[i]) printf("%04X ", codeword[i]);
		else printf("%04X ", data[i]);
	}
	printf("\n");

	for (int i=0; i<k; i++){//Check the correctness of the result
		if(erasure[i] == 1)
			if(data[i] != codeword[i]){
				printf("Decoding Error!\n");
				return;
			}
	}
	printf("Decoding is successful!\n");
	return;
}


#include <assert.h>

int flt_roundtrip() {
	const int N = 16;
	GFSymbol expected[16] = {
		1, 2, 3, 5, 8, 13, 21, 44,
		65, 0, 0xFFFF, 2, 3, 5, 7, 11,
	};
	GFSymbol data[N];
	memcpy(data, expected, N * sizeof(GFSymbol));


	FLT(data, N, N/4);
	printf("novel basis(c)\n");
	for(int i=0; i<N; i++){
		printf("%04X ", data[i]);
	}
	printf("\n");
	IFLT(data, N, N/4);
	for(int i=0; i<N; i++){
		assert(data[i] == expected[i]);
	}
	return 0;
}

int main(){
	// flt_roundtrip();

	init();//fill log table and exp table
	init_dec();//compute factors used in erasure decoder
	test(32);//test(int k), k: message size
	return 1;
}
