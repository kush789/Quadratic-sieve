#include <stdio.h>
#include "flint.h"
#include <math.h>

/* B = e ^ ((0.5)*((ln(n))*(ln(ln(n))))^0.5) */

mp_limb_t smoothness_bound(mp_limb_t n)
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


int main()
{
	return 0;
}