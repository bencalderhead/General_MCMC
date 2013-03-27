#ifndef GMCMC_BLOCKING_H
#define GMCMC_BLOCKING_H

#include <gmcmc/error.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * Type of parameter blocking.
 */
typedef enum {
  GMCMC_BLOCKING_TYPE_FIXED,    /*!< Fixed parameter blocking.  */
  GMCMC_BLOCKING_TYPE_RANDOM    /*!< Random parameter blocking. */
} gmcmc_blocking_type;

/*!
 * Parameter blocking.
 */
typedef struct __gmcmc_block_st gmcmc_block;

/*!
 * Creates a block.
 *
 * @param [out] block    the block to create
 * @param [in]  indices  the indices of parameters in the block
 * @param [in]  n        the number of parameters in the block
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_INVALID_ARGUMENT if n is negative,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to create
 *               the block.
 */
gmcmc_error gmcmc_block_create(gmcmc_block **, const int *, int);

/*!
 * Creates a block which is a deep copy of another block.
 *
 * @param [out] dest  the destination block
 * @param [in]  src   the source block
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to create
 *               the block.
 */
gmcmc_error gmcmc_block_create_copy(gmcmc_block **, const gmcmc_block *);

/*!
 * Performs a deep copy of one block into another block.
 *
 * @param [out] dest  the destination block
 * @param [in]  src   the source block
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_OUT_OF_MEMORY if there was not enough memory to create
 *               a copy of the contents of the block.
 */
gmcmc_error gmcmc_block_copy(gmcmc_block *, const gmcmc_block *);

/*!
 * Destroys a block.
 *
 * @param [in] block  the block to destroy.
 */
void gmcmc_block_destroy(gmcmc_block *);

/*!
 * Gets the number of parameters in a block.
 *
 * @param [in] block  the block
 *
 * @return the number of parameters in the block.
 */
int gmcmc_block_get_size(const gmcmc_block *);

/*!
 * Gets the index of one of the parameters in a block.
 *
 * @param [in]  block  the block
 * @param [in]  i      the index of the parameter in the block
 * @param [out] res    the index of the parameter in the model
 *
 * @return GMCMC_SUCCESS on success,
 *         GMCMC_ERROR_INVALID_ARGUMENT if the index is less than zero or greater
 *               than the size of the block.
 */
gmcmc_error gmcmc_block_get_index(const gmcmc_block *, int, int *);

#ifdef __cplusplus
}
#endif

#endif
