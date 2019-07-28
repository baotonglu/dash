#include "utils.h"
#include "libpmemobj.h"

static const char *pool_name = "pmem_hash.data";
static const char *layout_name = "hashtable";
static const uint32_t pool_size = 1024 * 1024 * 1024;

int main()
{
    PMEMobjpool *tmp_pool;
    if (!FileExists(pool_name))
    {
        tmp_pool = pmemobj_create(pool_name, layout_name, pool_size, CREATE_MODE_RW);
        if (tmp_pool == nullptr)
        {
            LOG_FATAL("failed to create a pool;");
        }
    }
    else
    {
        tmp_pool = pmemobj_open(pool_name, layout_name);
        if (tmp_pool == nullptr)
        {
            LOG_FATAL("failed to open the pool;");
        }
    }
}
