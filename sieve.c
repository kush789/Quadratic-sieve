#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include <math.h>

/* A table of prime numbers below 700 */

mp_limb_t prime_table[] = {   2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 
				 	  		 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 
					 		107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 
					 		167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 
					 		229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 
					 		283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 
					 		359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 
					 		431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 
					 		491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 
					 		571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 
					 		641, 643, 647, 653, 659, 661, 673, 677, 683, 691 };

					 /* the above table consists of 125 values */

/* B = e ^ ((0.5)*((ln(n))*(ln(ln(n))))^0.5) */

mp_limb_t
get_smoothness_bound(mp_limb_t n)
{
    mp_limb_t ret;
    double val, loge2, lnn, lnlnn;

    val = n;
    loge2 = 0.69314781;

    lnn = log2(val) * loge2;
    lnlnn = log2(lnn) * loge2;

    val = ceil(lnn * lnlnn);
    val = sqrt(val) * 0.5;
    val = pow(2.71, val);
    
    ret = ceil(val);
    return ret;
}

mp_limb_t*
get_primes(mp_limb_t B, mp_limb_t* num)
{
	int count, i;
	mp_limb_t* prime_arr;
	count = 0;
	i = 0;

	while(prime_table[i++]<=B)
		count+=1;

	prime_arr = (mp_limb_t*)malloc(count*sizeof(mp_limb_t));

	for (i = 0;i<count;i++)
		prime_arr[i] = prime_table[i];
	*num = count;
	return prime_arr;
}

mp_limb_t*
get_factor_arr(mp_limb_t n, mp_limb_t prime_arr[], mp_limb_t num, mp_limb_t* fac_count)
{
	mp_limb_t ninv;
	mp_limb_t* factor_base;
	int legendre_value[num], i, j, count;
	count = 0;
	j = 0;

	for (i = 0;i<num;i++)			/* calculating legendre symbol values */
	{
		ninv = n_preinvert_limb(prime_arr[i]);
		legendre_value[i] = n_powmod2_ui_preinv(n, (prime_arr[i] - 1)/2, prime_arr[i], ninv);
	}

	for (i = 0;i<num;i++)			/* calculating number of primes for which (a|p) = 1 */
		if (legendre_value[i] == 1)
			count+=1;

	factor_base = (mp_limb_t*)malloc(count*sizeof(mp_limb_t));

	for (i = 0;i<num;i++)			/* populating factor base */
		if (legendre_value[i] == 1)
			factor_base[j++] = prime_arr[i];

	*fac_count = count;
	return factor_base;
}

mp_limb_t*
get_factor_base(mp_limb_t n, mp_limb_t* factor_base_count)
{
	mp_limb_t B = get_smoothness_bound(n);
	mp_limb_t prime_count, factor_count;
	mp_limb_t* prime_arr = get_primes(B, &prime_count);
	mp_limb_t* factor_base = get_factor_arr(n, prime_arr, prime_count, &factor_count);
	*factor_base_count = factor_count;

	return factor_base;
}

int main()
{
	return 0;
}