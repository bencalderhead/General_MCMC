#include "mcmc/blocking.h"

/**
 * MCMC parameter block structure.
 */
struct __mcmc_block_st {
  int * indices;        // Indices of parameters in the block
  int n;                // Number of parameters in the block
};

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
mcmc_error mcmc_block_create(mcmc_block * block, const int * indices, int n) {
  if (n < 0)
    return MCMC_ERROR_INVALID_ARGUMENT;

  if ((*block = malloc(sizeof(struct __mcmc_block_st))) == NULL)
    return MCMC_ERROR_OUT_OF_MEMORY;

  if (((*block)->indices = malloc((size_t)n * sizeof(int))) == NULL) {
    free(*block);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  for (int i = 0; i < n; i++) {
    if (indices[i] < 0) {
      free((*block)->indices);
      free(*block);
      return MCMC_ERROR_INVALID_ARGUMENT;
    }
    (*block)->indices[i] = indices[i];
  }

  return MCMC_SUCCESS;
}

/**
 * Destroys a block.
 *
 * @param block  the block to destroy
 */
void mcmc_block_destroy(mcmc_block block) {
  free(block->indices);
  free(block);
}

int mcmc_block_get_size(const mcmc_block block) {
  return block->n;
}

int mcmc_block_get_index(const mcmc_block block, int i) {
  return block->indices[i];
}

/**
 * MCMC parameter blocking structure.
 */
struct __mcmc_blocking_st {
  mcmc_blocking_type type;      // Blocking type
  const mcmc_block * blocks;    // Blocks
  int n;                        // Number of blocks
};

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
mcmc_error mcmc_blocking_create(mcmc_blocking * blocking, mcmc_blocking_type type,
                                const mcmc_block * blocks, int n) {
  if (n < 0)
    return MCMC_ERROR_INVALID_ARGUMENT;

  if ((*blocking = malloc(sizeof(struct __mcmc_blocking_st))) == NULL)
    return MCMC_ERROR_OUT_OF_MEMORY;

  (*blocking)->type = type;

  if (((*blocking)->blocks = malloc((size_t)n * sizeof(mcmc_block))) == NULL) {
    free(*blocking);
    return MCMC_ERROR_OUT_OF_MEMORY;
  }

  for (int i = 0; i < n; i++)
    (*blocking)->blocks[i] = blocks[i];

  return MCMC_SUCCESS;
}

/**
 * Destroys an MCMC parameter blocking.
 *
 * @param blocking  the parameter blocking to destroy.
 */
void mcmc_blocking_destroy(mcmc_blocking blocking) {
  free(blocking->blocks);
  free(blocking);
}

mcmc_blocking_type mcmc_blocking_get_type(const mcmc_blocking blocking) {
  return blocking->type;
}

int mcmc_blocking_get_num_blocks(const mcmc_blocking blocking) {
  return blocking->n;
}

mcmc_block mcmc_blocking_get_block(const mcmc_blocking blocking, int i) {
  return blocking->blocks[i];
}
