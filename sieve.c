#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
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
    
    ret = ceil(val) + 20;
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

	prime_arr = (mp_limb_t*)flint_malloc(count*sizeof(mp_limb_t));

	for (i = 0;i<count;i++)
		prime_arr[i] = prime_table[i];
	*num = count;
	return prime_arr;
}

fmpz*
get_factor_arr(mp_limb_t n, mp_limb_t prime_arr[], mp_limb_t num, mp_limb_t* fac_count)
{
	mp_limb_t ninv;
	fmpz* factor_base;
	int legendre_value[num], i, j, count;
	count = 0;
	j = 1;

	for (i = 0;i<num;i++)			/* calculating legendre symbol values */
	{
		ninv = n_preinvert_limb(prime_arr[i]);
		legendre_value[i] = n_powmod2_ui_preinv(n, (prime_arr[i] - 1)/2, prime_arr[i], ninv);
	}

	for (i = 0;i<num;i++)			/* calculating number of primes for which (a|p) = 1 */
		if (legendre_value[i] == 1)
			count+=1;
	count+=1;
	factor_base = (fmpz*)flint_malloc(count*sizeof(mp_limb_t));
	factor_base[0] = -1;
	for (i = 0;i<num;i++)			/* populating factor base */
		if (legendre_value[i] == 1)
			factor_base[j++] = prime_arr[i];

	*fac_count = count;
	return factor_base;
}

fmpz*
get_factor_base(mp_limb_t n, mp_limb_t* factor_base_count)
{
	mp_limb_t B = get_smoothness_bound(n);
	mp_limb_t prime_count, factor_count;
	mp_limb_t* prime_arr = get_primes(B, &prime_count);
	fmpz* factor_base = get_factor_arr(n, prime_arr, prime_count, &factor_count);
	*factor_base_count = factor_count;

	return factor_base;
}

void
get_sieveing_range(mp_limb_t n, mp_limb_t* min, mp_limb_t* max)
{
	mp_limb_t B, M, sqrtn;

	sqrtn = n_sqrt(n);
	B = get_smoothness_bound(n);

	*min = sqrtn - B;
	*max = sqrtn + B;
}

int**
get_relation_matrix(mp_limb_t n, fmpz* factor_base, mp_limb_t factor_base_count, mp_limb_t* numbers, mp_limb_t* numbers_count)
{
	int i, j, count, B_smooth_count;
	mp_limb_t min, max, range, iter, currdiff;
	fmpz* diff_values;

	fmpz_t x, xsq, currfactor;
	fmpz_init(x);
	fmpz_init(xsq);
	fmpz_init(currfactor);
	B_smooth_count = 0;

	get_sieveing_range(n, &min, &max);
	range = max - min + 1;

	/*  table containing x^2 - n values for x in range [min, max] */

	diff_values = (fmpz*)flint_malloc(range*sizeof(mp_limb_t));
	for (iter = min; iter<=max;iter+=1)
	{
		fmpz_set_ui(x, iter);
		*x = (*x)*(*x);
		*x = *x - n;
		diff_values[iter-min] = *x;
	}

	int factorization[range][factor_base_count];

	for (i = 0;i<range;i++)						/* checking if negative, for factor -1 */
		if (diff_values[i] < 0)
		{
			factorization[i][0] = 1;
			diff_values[i]*=(-1);
		}
		else
			factorization[i][0] = 0;

	/* Checking for B smoothness for positive factors */

	for (i = factor_base_count -1;i>0;i--)		/* for each prime in factor base */
	{
		*currfactor = factor_base[i];

		for (j = 0;j<range;j++)					/* for each element in array diff_values */
		{
			count = 0;
				while (1)
				{
					if (diff_values[j] == 1)
						break;					/* if factored completely */
						
					if (diff_values[j]%(*currfactor) == 0)
					{
						count+=1;				/* count powers of each prime factor */
						diff_values[j]/=(*currfactor);
					}
					else
						break;
				}

			factorization[j][i] = count%2;		
		}
	}

	for (j = 0;j<range;j++)						/* Counting numbers which are B smooth */
		if (diff_values[j] == 1)
			B_smooth_count+=1;
	
	/* array to store B smooth numbers */
		
	mp_limb_t* B_smooth_numbers = (mp_limb_t*)flint_malloc(B_smooth_count * sizeof(mp_limb_t));
	j = 0;
	for (i = 0; i < range;i++)
	{
		if (diff_values[i] == 1)
			B_smooth_numbers[j++] = min+i;		/* storing B smooth numbers */
	}

	int** sieve_matrix;							/* declaring space for sieve matrix */
	sieve_matrix = (int**) flint_malloc(factor_base_count * sizeof(int*));
	for (i = 0;i<factor_base_count;i++)
		sieve_matrix[i] = (int*) flint_malloc(B_smooth_count * sizeof(int));

	for (i = 0;i<factor_base_count;i++)			/* row ->factors , columns -> B smooth numbers */
		for (j = 0;j<B_smooth_count;j++)
			sieve_matrix[i][j] = factorization[B_smooth_numbers[j]-min][i];


	numbers = B_smooth_numbers;
	*numbers_count = B_smooth_count;
	return sieve_matrix;
}

int main()
{
	return 0;
}