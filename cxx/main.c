
#include "./RSErasureCode.h"
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
	return roundtrip(32, 4);//test(n, k), k: message size, n: domain size
}
