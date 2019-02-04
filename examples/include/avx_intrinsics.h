#ifndef _AVX_INTRIN_RACE_
#define _AVX_INTRIN_RACE_

#if (VECLEN == 4)
    #define VECREG __m256d
    #define SET_ZERO _mm256_setzero_pd()
    #define BROADCAST(addr) _mm256_broadcast_sd(addr)
    #define BROADCAST_128(addr) _mm256_broadcast_pd(addr)
    #define CAST_TO_128(val) _mm256_castpd256_pd128(val)
    #define SET_FROM_128(val1, val2) _mm256_set_m128d(val1, val2)

    #define LOAD(addr) _mm256_loadu_pd(addr)
    #define STORE(addr,val) _mm256_storeu_pd(addr,val)
    #define ADD(val1, val2) _mm256_add_pd(val1,val2)
    #define MUL(val1, val2) _mm256_mul_pd(val1, val2)
    #define FMA(val1, val2, val3) ADD(MUL(val1, val2),val3)

    #define COMPLEX_MUL(a,b) _mm256_addsub_pd(_mm256_mul_pd(_mm256_shuffle_pd(b,b,0),a),_mm256_mul_pd(_mm256_shuffle_pd(b,b,0xF),_mm256_shuffle_pd(a,a,5)))

    #define COMPLEX_MUL_CONJ(b,a) _mm256_addsub_pd(_mm256_mul_pd(_mm256_shuffle_pd(b,b,0),a),_mm256_mul_pd(_mm256_mul_pd(_mm256_shuffle_pd(b,b,0xF),_mm256_set1_pd(-1.)),_mm256_shuffle_pd(a,a,5)))

    //val_1*val_2 + val_3
#elif (VECLEN == 8)
    #define VECREG __m512d
    #define SET_ZERO _mm512_setzero_pd()
    #define BROADCAST(addr) _mm512_set1_pd(*addr)
    #define BROADCAST_128(addr) _mm512_broadcast_pd(*addr)
    #define CAST_TO_128(val) _mm512_castpd512_pd128(val)
    #define SET_FROM_128(val1, val2) _mm512_broadcast_f64x4(_mm256_set_m128d(val1,val2))
    #define LOAD(addr) _mm512_loadu_pd(addr)
    #define STORE(addr, val) _mm512_storeu_pd(addr, val)
    #define ADD(val1, val2) _mm512_add_pd(val1, val2)
    #define MUL(val1, val2) _mm512_mul_pd(val1, val2)
    #define FMA(val1, val2, val3) _mm512_fmadd_pd(val1, val2, val3)

    #define COMPLEX_MUL(a, b) _mm512_fmaddsub_pd(_mm512_shuffle_pd(b,b,0),a, _mm512_mul_pd(_mm512_shuffle_pd(b,b,0xF), _mm512_shuffle_pd(a,a,85)))

    #define COMPLEX_MUL_CONJ(b, a) _mm512_fmaddsub_pd(_mm512_shuffle_pd(b,b,0),a,(_mm512_fnmadd_pd(_mm512_shuffle_pd(b,b,255),_mm512_shuffle_pd(a,a,85), _mm512_setzero_pd())))

#endif

#endif
