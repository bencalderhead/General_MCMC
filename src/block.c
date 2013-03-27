#include <gmcmc/block.h>
#include <stdlib.h>
#include <string.h>

/*!
 * Parameter blocking.
 */
struct __gmcmc_block_st {
  int * indices;        /**< Indices of parameters in the block */
  int n;                /**< Number of parameters in the block  */
};

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
gmcmc_error gmcmc_block_create(gmcmc_block ** block, const int * indices, int n) {
  if (n <= 0)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);

  // Allocate space for the block
  if ((*block = malloc(sizeof(gmcmc_block))) == NULL)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);

  // Allocate space for the block indices
  if (((*block)->indices = malloc((size_t)n * sizeof(int))) == NULL) {
    free(*block);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
  }

  // Copy the indices into place
  (*block)->indices = memcpy((*block)->indices, indices, (size_t)n * sizeof(int));
  (*block)->n = n;

  return GMCMC_SUCCESS;
}

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
gmcmc_error gmcmc_block_create_copy(gmcmc_block ** dest, const gmcmc_block * src) {
  // Allocate space for the block
  if ((*dest = malloc(sizeof(gmcmc_block))) == NULL)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);

  // Allocate space for the block indices
  if (((*dest)->indices = malloc((size_t)src->n * sizeof(int))) == NULL) {
    free(*dest);
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);
  }

  // Copy the block indices from the source block to the destination
  (*dest)->indices = memcpy((*dest)->indices, src->indices, (size_t)src->n * sizeof(int));
  (*dest)->n = src->n;

  return GMCMC_SUCCESS;
}

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
gmcmc_error gmcmc_block_copy(gmcmc_block * dest, const gmcmc_block * src) {
  // If dest == src return now
  if (dest == src)
    return GMCMC_SUCCESS;

  // Resize the block only if needed (possibility of failure)
  if (dest->n != src->n) {
    // Allocate new space for the indices so that the original block is untouched on failure
    int * indices;
    if ((indices = malloc((size_t)src->n * sizeof(int))) == NULL)
      GMCMC_ERROR_HANDLER(GMCMC_ERROR_OUT_OF_MEMORY);

    // Update the destination block
    free(dest->indices);
    dest->indices = indices;
    dest->n = src->n;
  }

  // Copy the block indices from the source block to the destination
  dest->indices = memcpy(dest->indices, src->indices, (size_t)src->n * sizeof(int));

  return GMCMC_SUCCESS;
}

/*!
 * Destroys a block.
 *
 * @param [in] block  the block to destroy.
 */
void gmcmc_block_destroy(gmcmc_block * block) {
  free(block->indices);
  free(block);
}

/*!
 * Gets the number of parameters in a block.
 *
 * @param [in] block  the block
 *
 * @return the number of parameters in the block.
 */
int gmcmc_block_get_size(const gmcmc_block * block) {
  return block->n;
}

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
gmcmc_error gmcmc_block_get_index(const gmcmc_block * block, int i, int * res) {
  if (i < 0 || i >= block->n)
    GMCMC_ERROR_HANDLER(GMCMC_ERROR_INVALID_ARGUMENT);
  *res = block->indices[i];
  return GMCMC_SUCCESS;
}
