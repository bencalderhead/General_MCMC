#ifndef MEMPOOL_H
#define MEMPOOL_H

#include <stddef.h>

/**
 * A block of memory.
 */
typedef struct __block_st block;

/**
 * Gets an aligned pointer into the memory allocated in the block.
 *
 * @param b  the block.
 *
 * @return an aligned pointer.
 */
void * block_ptr(block *);

/**
 * Gets an aligned constant pointer into the memory allocated in the block.
 *
 * @param b  the block.
 *
 * @return an aligned constant pointer.
 */
const void * block_ptr_const(const block *);


/**
 * A pool of fixed-size memory blocks.
 */
typedef struct __mempool_st mempool;

/**
 * Allocates a pool of memory blocks.
 *
 * @param n  the size of the blocks to allocate.
 *
 * @return a pointer to the newly allocated memory pool, or NULL if there is not
 * enough memory to create the pool.
 */
mempool * mempool_alloc(size_t);

/**
 * Frees all memory blocks managed by the pool, then frees the pool.
 *
 * @param pool  the pool to free.
 */
void mempool_free(mempool *);

/**
 *  Allocates a block from the pool.
 *
 * @param pool  the pool to allocate from.
 *
 * @return the allocated block, or NULL if a block could not be allocated.
 */
block * block_alloc(mempool *);

/**
 * Returns a block to the pool so that it may be reallocated.
 *
 * @param b  the block to return.
 */
void block_free(block *);

#endif

//#include "mempool.h"
#include <stdlib.h>

// Alignment requirement in bytes
static const size_t alignment = 16u;

struct __block_st {
  void * base, * aligned;       // Base pointer, aligned pointer
  block * prev, * next;         // Pointer to previous and next blocks
  mempool * pool;               // The pool this block belongs to
};

void * block_ptr(block * b) { return b->aligned; }
const void * block_ptr_const(const block * b) { return b->aligned; }

typedef struct __mempool_st {
  block * allocated;
  block * available;
  size_t size;
} mempool;

mempool * mempool_alloc(size_t n) {
  mempool * pool = malloc(sizeof(struct __mempool_st));

  if (pool != NULL) {
    pool->allocated = NULL;
    pool->available = NULL;
    pool->size = n;
  }

  return pool;
}

void mempool_free(mempool * pool) {
  // Free the allocated blocks
  while (pool->allocated != NULL) {
    block * b = pool->allocated;
    pool->allocated = b->next;
    free(b->base);
    free(b);
  }

  // Free the available blocks
  while (pool->available != NULL) {
    block * b = pool->available;
    pool->available = b->next;
    free(b->base);
    free(b);
  }

  free(pool);
}

block * block_alloc(mempool * pool) {
  block * b;

  // If there are blocks available
  if (pool->available != NULL) {
    // Remove the block from the list of available blocks
    b = pool->available;
    pool->available = pool->available->next;
  }
  else {
    // Allocate a new block from the system
    b = malloc(sizeof(block));
    if (b == NULL)
      return NULL;

    b->base = malloc(pool->size + alignment);
    if (b->base == NULL) {
      free(b);
      return NULL;
    }

    b->aligned = (void *)((((size_t)b->base) + alignment - 1u) & ~(alignment - 1u));
    b->pool = pool;
  }

  // Add the block to the list of blocks in use
  b->next = pool->allocated;
  pool->allocated = b;
  b->prev = NULL;
  if (b->next != NULL)
    b->next->prev = b;

  return b;
}

void block_free(block * b) {
  // Remove the block from the list of allocated blocks
  if (b->pool->allocated == b)
    b->pool->allocated = b->next;
  else
    b->prev->next = b->next;
  if (b->next != NULL)
    b->next->prev = b->prev;

  // Add the block to list of available blocks
  b->next = b->pool->available;
  b->pool->available = b;
  b->prev = NULL;
  if (b->next != NULL)
    b->next->prev = b;
}

int main() {
  mempool * pool;
  block * a, * b, * c, * d;

  pool = mempool_alloc(sizeof(double));

  a = block_alloc(pool);
  b = block_alloc(pool);
  c = block_alloc(pool);
  block_free(c);
  d = block_alloc(pool);
  c = block_alloc(pool);
  block_free(a);
  block_free(b);
  block_free(c);

  mempool_free(pool);

  return 0;
}
