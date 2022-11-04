/*
 * Write a parallel program that computes sums to arbitrary
 * precision after the decimal point
 * Use GNU 
 * ./hps_gmp <n> <d>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "mpi.h"

#define PREC 256L // set precision 

#define NLIMBS(n) \
    ((mp_size_t)((__GMP_MAX(53, n) + 2 * GMP_NUMB_BITS - 1) / GMP_NUMB_BITS))
#define MPF_D_SIZE(prec) (((prec) + 1) * sizeof(mp_limb_t))
#define MPZ_D_SIZE(alloc) ((alloc) * sizeof(mp_limb_t))
#define MPF_PACKED_LIMBS(prec) ((prec) + 4);
#define MPZ_PACKED_LIMBS(alloc) ((alloc) + 2);
#define MPF_PACKED_BYTES ((NLIMBS(mpf_get_default_prec()) + 3) * sizeof(mp_limb_t))

mp_limb_t *mpf_pack(mp_limb_t *dest, mpf_t *src, int n)
{
    if (dest == NULL)
    {
        size_t size = 0;
        for (int i = 0; i < n; i++)
            size += MPF_PACKED_LIMBS(src[i]->_mp_prec);
        dest = (mp_limb_t *)malloc(size * sizeof(mp_limb_t));
    }
    size_t offset = 0;
    for (int i = 0; i < n; i++)
    {
        int prec = dest[offset] = src[i]->_mp_prec;
        dest[offset + 1] = src[i]->_mp_size;
        dest[offset + 2] = src[i]->_mp_exp;
        memcpy(&dest[offset + 3], src[i]->_mp_d, MPF_D_SIZE(prec));
        offset += MPF_PACKED_LIMBS(prec);
    }
    return dest;
}

mpf_t *mpf_unpack(mpf_t *dest, mp_limb_t *src, int n)
{
    if (dest == NULL)
    {
        dest = (mpf_t *)malloc(sizeof(mpf_t) * n);
    }
    else
    {
        for (int i = 0; i < n; i++)
            if (dest[i]->_mp_d)
                free(dest[i]->_mp_d);
    }
    size_t offset = 0;
    for (int i = 0; i < n; i++)
    {
        int prec = dest[i]->_mp_prec = src[offset];
        dest[i]->_mp_size = src[offset + 1];
        dest[i]->_mp_exp = src[offset + 2];
        dest[i]->_mp_d = malloc(MPF_D_SIZE(prec));
        memcpy(dest[i]->_mp_d, &src[offset + 3], MPF_D_SIZE(prec));
        offset += MPF_PACKED_LIMBS(prec);
    }
    return dest;
}

mp_limb_t *mpz_pack(mp_limb_t *dest, mpz_t *src, int n)
{
    if (dest == NULL)
    {
        size_t size = 0;
        for (int i = 0; i < n; i++)
            size += MPZ_PACKED_LIMBS(src[i]->_mp_alloc);
        dest = (mp_limb_t *)malloc(size * sizeof(mp_limb_t));
    }
    size_t offset = 0;
    for (int i = 0; i < n; i++)
    {
        int alloc = dest[offset] = src[i]->_mp_alloc;
        dest[offset + 1] = src[i]->_mp_size;
        memcpy(&dest[offset + 2], src[i]->_mp_d, MPZ_D_SIZE(alloc));
        offset += MPZ_PACKED_LIMBS(alloc);
    }
    return dest;
}

mpz_t *mpz_unpack(mpz_t *dest, mp_limb_t *src, int n)
{
    if (dest == NULL)
    {
        dest = (mpz_t *)malloc(sizeof(mpz_t) * n);
    }
    else
    {
        for (int i = 0; i < n; i++)
            if (dest[i]->_mp_d)
                free(dest[i]->_mp_d);
    }
    size_t offset = 0;
    for (int i = 0; i < n; i++)
    {
        int alloc = dest[i]->_mp_alloc = src[offset];
        dest[i]->_mp_size = src[offset + 1];
        dest[i]->_mp_d = malloc(MPZ_D_SIZE(alloc));
        memcpy(dest[i]->_mp_d, &src[offset + 2], MPZ_D_SIZE(alloc));
        offset += MPZ_PACKED_LIMBS(alloc);
    }
    return dest;
}

void mpf_packed_add(void *_in, void *_inout, int *len, MPI_Datatype *datatype)
{
    mpf_t *in = mpf_unpack(NULL, (mp_limb_t *)_in, *len);
    mpf_t *inout = mpf_unpack(NULL, (mp_limb_t *)_inout, *len);

    for (int i = 0; i < *len; i++)
        mpf_add(inout[i], in[i], inout[i]);

    mpf_pack((mp_limb_t *)_inout, inout, *len);

    for (int i = 0; i < *len; i++)
    {
        mpf_clear(in[i]);
        mpf_clear(inout[i]);
    }
    free(in);
    free(inout);
}

void mpz_packed_add(void *_in, void *_inout, int *len, MPI_Datatype *datatype)
{
    mpz_t *in = mpz_unpack(NULL, (mp_limb_t *)_in, *len);
    mpz_t *inout = mpz_unpack(NULL, (mp_limb_t *)_inout, *len);

    for (int i = 0; i < *len; i++)
        mpz_add(inout[i], in[i], inout[i]);

    mpz_pack((mp_limb_t *)_inout, inout, *len);

    for (int i = 0; i < *len; i++)
    {
        mpz_clear(in[i]);
        mpz_clear(inout[i]);
    }
    free(in);
    free(inout);
}

int main(int argc, char *argv[])
{
    int id;    /* Process rank */
    int p;     /* Number of processes */
    int n = 1; // summation n terms
    int d = 1; // digits of decimal

    mpf_t sum;
    mpf_t global_sum;
    mpf_t one;
    mpf_t tempf;
    double elapsed_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    mpf_set_default_prec(PREC);
    // Create an MPI Datatype based on the target precision
    MPI_Datatype MPI_MPF;
    MPI_Type_contiguous(MPF_PACKED_BYTES, MPI_CHAR, &MPI_MPF);
    MPI_Type_commit(&MPI_MPF);

    // Create an MPI Operation for adding mpfs
    MPI_Op MPI_SUM_MPF;
    MPI_Op_create(mpf_packed_add, 1, &MPI_SUM_MPF);

    // Create the number 1

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    n = atoi(argv[1]);
    d = atoi(argv[2]);

    // mpf_set_default_prec(65536UL);
    mpf_init(sum);
    mpf_init(global_sum);
    mpf_init(one);
    mpf_init(tempf);

    mpf_set_ui(one, 1UL);
    mpf_set_ui(sum, 0UL);

    // boardcast n, d to process
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // printf("n = %d, d = %d\n", n, d);

    for (int i = id + 1; i <= n; i += p)
    {
        mpf_div_ui(tempf, one, i);
        mpf_add(sum, sum, tempf);
    }
    mpf_clear(one);
    mpf_clear(tempf);

    mp_limb_t *packed = mpf_pack(NULL, &sum, 1);
    mp_limb_t *sum_packed = malloc(MPF_PACKED_BYTES);
    MPI_Reduce(packed, sum_packed, 1, MPI_MPF, MPI_SUM_MPF, 0, MPI_COMM_WORLD);
    elapsed_time += MPI_Wtime();


    if (id == 0)
    {
        mpf_unpack(&global_sum, sum_packed, 1);
        mpf_out_str(stdout, 10, d, global_sum);
        printf("\nElapsed time : %lf\n", elapsed_time);
    }
    mpf_clear(sum);
    mpf_clear(global_sum);
    free(packed);
    free(sum_packed);

    MPI_Finalize();
    return 0;
}
