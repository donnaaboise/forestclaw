#include "fc2d_multigrid_vector.hpp"
#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>
#include <fclaw2d_patch.h>
#include "fc2d_multigrid.h"
#include "fc2d_multigrid_options.h"

using namespace ThunderEgg;

/**
 * @brief Get the number of local cells in the forestclaw domain
 * 
 * @param glob the forestclaw glob
 * @return int the number of local cells
 */
static int get_num_local_cells(fclaw2d_global_t* glob){
    fclaw2d_clawpatch_options_t *clawpatch_opt = fclaw2d_clawpatch_get_options(glob);

    fclaw2d_domain_t *domain = glob->domain;

    return domain->local_num_patches * clawpatch_opt->mx * clawpatch_opt->my;
}

fc2d_multigrid_vector::fc2d_multigrid_vector(fclaw2d_global_t *glob)
    : Vector<2>( MPI_COMM_WORLD,
                 fclaw2d_clawpatch_get_options(glob)->meqn,
                 glob->domain->local_num_patches, 
                 get_num_local_cells(glob) ) 
{

    fclaw2d_clawpatch_options_t *clawpatch_opt =
        fclaw2d_clawpatch_get_options(glob);

    fclaw2d_domain_t *domain = glob->domain;

    meqn = clawpatch_opt->meqn;

    ns[0] = clawpatch_opt->mx;
    ns[1] = clawpatch_opt->my;

    mbc = clawpatch_opt->mbc;

    strides[0] = 1;
    strides[1] = ns[0] + 2 * mbc;
    eqn_stride = strides[1] * (ns[1] + 2 * mbc);
    /*
    eqn_stride = 1;
    strides[0] = clawpatch_opt->meqn;
    strides[1] = strides[0] * (ns[0] + 2 * mbc);
    */

    patch_data.resize(domain->local_num_patches * clawpatch_opt->meqn);
    fclaw2d_global_iterate_patches(glob, enumeratePatchData, this);
}
void fc2d_multigrid_vector::enumeratePatchData(fclaw2d_domain_t *domain,
                                               fclaw2d_patch_t *patch,
                                               int blockno, int patchno,
                                               void *user) 
{
    fclaw2d_global_iterate_t *g = (fclaw2d_global_iterate_t *)user;

    int meqn;

    int global_num, local_num, level;

    fc2d_multigrid_vector &vec = *(fc2d_multigrid_vector *)g->user;

    fclaw2d_patch_get_info(domain, patch, blockno, patchno, &global_num,
                           &local_num, &level);

    fclaw2d_clawpatch_soln_data(g->glob, patch, &vec.patch_data[local_num],
                                &meqn);
    // update point to point to non ghost cell
    for(int c = 0; c < meqn; c++){
        vec.patch_data[local_num * meqn + c] += vec.eqn_stride * c +
                                            vec.mbc * vec.strides[0] +
                                            vec.mbc * vec.strides[1];
    }
}
LocalData<2> fc2d_multigrid_vector::getLocalDataPriv(int component_index, int local_patch_id) const 
{
    return LocalData<2>(patch_data[local_patch_id * meqn + component_index], strides, ns,mbc);
}
LocalData<2> fc2d_multigrid_vector::getLocalData(int component_index, int local_patch_id) 
{
    return getLocalDataPriv(component_index, local_patch_id);
}
const LocalData<2> fc2d_multigrid_vector::getLocalData(int component_index, int local_patch_id) const 
{
    return getLocalDataPriv(component_index, local_patch_id);
}
