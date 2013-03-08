#ifndef MCMC_BLOCKING_H
#define MCMC_BLOCKING_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * MCMC parameter block structure.
 */
typedef struct __mcmc_block_st * mcmc_block;

/**
 * Creates a new block.
 *
 * @param block    the block to create
 * @param indices  indices of parameters in the block
 * @param n        the number of parameters in the block
 *
 * @return MCMC_SUCCESS on success,
 *         MCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create
 *                    another block,
 *         MCMC_ERROR_INVALID_ARGUMENT if any of the elements of indices are
 *                    negative or if n is negative.
 */
mcmc_error mcmc_block_create(mcmc_block *, int *, int);

/**
 * Creates a new block which is a copy of another block.
 *
 * @param dest  the block to create
 * @param src   the block to copy from
 *
 * @return MCMC_SUCCESS on success,
 *         MCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create
 *                    another block.
 */
mcmc_error mcmc_block_clone(mcmc_block *, const mcmc_block);

/**
 * Copies one block into another
 *
 * @param dest  the block to copy into
 * @param src   the block to copy from
 *
 * @return the destination block.
 */
mcmc_block mcmc_block_copy(mcmc_block, const mcmc_block);

/**
 * Destroys a block.
 *
 * @param block  the block to destroy
 */
void mcmc_block_destroy(mcmc_block);

int mcmc_block_get_size(const mcmc_block);
int mcmc_block_get_index(const mcmc_block, int);


/**
 * Parameter blocking type.
 */
typedef enum { BLOCKING_FIXED, BLOCKING_RANDOM } mcmc_blocking_type;


/**
 * MCMC parameter blocking structure.
 */
typedef struct __mcmc_blocking_st * mcmc_blocking;

/**
 * Creates a new parameter blocking object.
 *
 * @param blocking  the blocking object to create
 * @param type      the type of blocking this object represents
 * @param blocks    the blocks in this parameter blocking
 * @param n         the number of blocks
 *
 * @return MCMC_SUCCESS on success,
 *         MCMC_ERROR_OUT_OF_MEMORY if there is not enough memory to create
 *                    another blocking object,
 *         MCMC_ERROR_INVALID_ARGUMENT if n is negative.
 */
mcmc_error mcmc_blocking_create(mcmc_blocking *, mcmc_blocking_type, const mcmc_block *, int);

/**
 * Destroys an MCMC parameter blocking.
 *
 * @param blocking  the parameter blocking to destroy.
 */
void mcmc_blocking_destroy(mcmc_blocking);

mcmc_blocking_type mcmc_blocking_get_type(const mcmc_blocking);
int mcmc_blocking_get_num_blocks(const mcmc_blocking);
mcmc_block mcmc_blocking_get_block(const mcmc_blocking, int);

#ifdef __cplusplus
}
#endif

#endif
