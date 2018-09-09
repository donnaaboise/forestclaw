/*
Copyright (c) 2012 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "fc2d_cudaclaw5.h"
#include "fc2d_cudaclaw5_fort.h"
#include "fc2d_cudaclaw5_options.h"


#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch.hpp>

#include <fclaw2d_clawpatch_output_ascii.h>
#include <fclaw2d_clawpatch_output_vtk.h>


#include <fclaw2d_patch.h>
#include <fclaw2d_global.h>
#include <fclaw2d_vtable.h>
#include <fclaw2d_defs.h>


static fc2d_cudaclaw5_vtable_t s_cudaclaw5_vt;

/* -------------------------- Clawpack solver functions ------------------------------ */

static
void cudaclaw5_setprob(fclaw2d_global_t *glob)
{
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();
    if (cuclaw5_vt->fort_setprob != NULL)
    {
        cuclaw5_vt->fort_setprob();
    }
}

/* This should only be called when a new fclaw2d_clawpatch_t is created. */
static
void cudaclaw5_setaux(fclaw2d_global_t *glob,
                      fclaw2d_patch_t *this_patch,
                      int this_block_idx,
                      int this_patch_idx)
{
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();
    if (cuclaw5_vt->fort_setaux == NULL)
    {
        return;
    }

    if (fclaw2d_patch_is_ghost(this_patch))
    {
        /* This is going to be removed at some point */
        return;
    }


    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux;


    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CUDACLAW5_SET_BLOCK(&this_block_idx);
    cuclaw5_vt->fort_setaux(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                          &maux,aux);
    CUDACLAW5_UNSET_BLOCK();
}

static
void cudaclaw5_qinit(fclaw2d_global_t *glob,
                     fclaw2d_patch_t *this_patch,
                     int this_block_idx,
                     int this_patch_idx)
{
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();
    
    FCLAW_ASSERT(cuclaw5_vt->fort_qinit != NULL); /* Must initialized */
    int mx,my,mbc,meqn,maux;
    double dx,dy,xlower,ylower;
    double *q, *aux;

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    CUDACLAW5_SET_BLOCK(&this_block_idx);
    cuclaw5_vt->fort_qinit(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
                         &maux,aux);
    CUDACLAW5_UNSET_BLOCK();
}

static
void cudaclaw5_b4step2(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t, double dt)

{
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();
    
    if (cuclaw5_vt->fort_b4step2 == NULL)
    {
        return;
    }

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CUDACLAW5_SET_BLOCK(&this_block_idx);
    cuclaw5_vt->fort_b4step2(&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                           &dx,&dy,&t,&dt,&maux,aux);
    CUDACLAW5_UNSET_BLOCK();
}

static
void cudaclaw5_src2(fclaw2d_global_t *glob,
                    fclaw2d_patch_t *this_patch,
                    int this_block_idx,
                    int this_patch_idx,
                    double t,
                    double dt)
{
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();
    
    if (cuclaw5_vt->fort_src2 == NULL)
    {
        return;
    }

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&q,&meqn);
    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    CUDACLAW5_SET_BLOCK(&this_block_idx);
    cuclaw5_vt->fort_src2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                        &dx,&dy,q,&maux,aux,&t,&dt);
    CUDACLAW5_UNSET_BLOCK();
}

static
void cudaclaw5_bc2(fclaw2d_global_t *glob,
                   fclaw2d_patch_t *this_patch,
                   int this_block_idx,
                   int this_patch_idx,
                   double t,
                   double dt,
                   int intersects_phys_bdry[],
                   int time_interp)
{
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();

    fc2d_cudaclaw5_options_t *clawopt = fc2d_cudaclaw5_get_options(glob);
    
    FCLAW_ASSERT(cuclaw5_vt->fort_bc2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    int *block_mthbc = clawopt->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[4];
    for(int i = 0; i < 4; i++)
    {
        if (intersects_phys_bdry[i])
        {
            mthbc[i] = block_mthbc[i];
        }
        else
        {
            mthbc[i] = -1;
        }
    }

    /*
      We may be imposing boundary conditions on time-interpolated data;
      and is being done just to help with fine grid interpolation.
      In this case, this boundary condition won't be used to update
      anything
    */
    fclaw2d_clawpatch_timesync_data(glob,this_patch,time_interp,&q,&meqn);

    CUDACLAW5_SET_BLOCK(&this_block_idx);
    cuclaw5_vt->fort_bc2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                   &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    CUDACLAW5_UNSET_BLOCK();

}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
static
double cudaclaw5_step2(fclaw2d_global_t *glob,
                       fclaw2d_patch_t *this_patch,
                       int this_block_idx,
                       int this_patch_idx,
                       double t,
                       double dt)
{
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();

    fc2d_cudaclaw5_options_t* cudaclaw_options;
    int level;
    double *qold, *aux;
    int mx, my, meqn, maux, mbc;
    double xlower, ylower, dx,dy;

    FCLAW_ASSERT(cuclaw5_vt->fort_rpn2 != NULL);
    FCLAW_ASSERT(cuclaw5_vt->fort_rpt2 != NULL);

    cudaclaw_options = fc2d_cudaclaw5_get_options(glob);

    // cp = fclaw2d_clawpatch_get_cp(this_patch);

    level = this_patch->level;

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);

    fclaw2d_clawpatch_save_current_step(glob, this_patch);

    fclaw2d_clawpatch_grid_data(glob,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(glob,this_patch,&qold,&meqn);

    int mwaves = cudaclaw_options->mwaves;

    int maxm = fmax(mx,my);

    double cflgrid = 0.0;

    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];


    int ierror = 0;
    //fc2d_cudaclaw5_flux2_t flux2 = cudaclaw_options->use_fwaves ?
    //                                CUDACLAW5_FLUX2FW : CUDACLAW5_FLUX2;
    cudaclaw5_fort_flux2_t flux2 = CUDACLAW5_FLUX2;
    int* block_corner_count = fclaw2d_patch_block_corner_count(glob,this_patch);
    CUDACLAW5_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, cudaclaw_options->method,
                          cudaclaw_options->mthlim, &cudaclaw_options->mcapa,
                          &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
                          work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
                          cuclaw5_vt->fort_rpn2, cuclaw5_vt->fort_rpt2,flux2,
                          block_corner_count, &ierror);

    FCLAW_ASSERT(ierror == 0);

    cudaclaw5_test();

    //CUDACLAW5_STEP2(&maxm, &meqn, &maux, &mbc, &mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
    //                fm, fp, gm, gp, cuclaw5_vt->rpn2, cuclaw5_vt->rpt2);
    /* Accumulate fluxes needed for conservative fix-up */
    if (cuclaw5_vt->fort_fluxfun != NULL)
    {
        /* Accumulate fluxes */
    }


    delete [] fp;
    delete [] fm;
    delete [] gp;
    delete [] gm;

    delete [] work;

    return cflgrid;
}

static
double cudaclaw5_update(fclaw2d_global_t *glob,
                        fclaw2d_patch_t *this_patch,
                        int this_block_idx,
                        int this_patch_idx,
                        double t,
                        double dt)
{    
    fc2d_cudaclaw5_vtable_t*  cuclaw5_vt = fc2d_cudaclaw5_vt();

    const fc2d_cudaclaw5_options_t* cudaclaw_options;
    cudaclaw_options = fc2d_cudaclaw5_get_options(glob);

    if (cuclaw5_vt->b4step2 != NULL)
    {
        cuclaw5_vt->b4step2(glob,
                          this_patch,
                          this_block_idx,
                          this_patch_idx,t,dt);
    }
    double maxcfl = cudaclaw5_step2(glob,
                                    this_patch,
                                    this_block_idx,
                                    this_patch_idx,t,dt);

    if (cudaclaw_options->src_term > 0 && cuclaw5_vt->src2 != NULL)
    {
        cuclaw5_vt->src2(glob,
                       this_patch,
                       this_block_idx,
                       this_patch_idx,t,dt);
    }
    return maxcfl;
}

/* ---------------------------------- Output functions -------------------------------- */

static
void cudaclaw5_output(fclaw2d_global_t *glob, int iframe)
{
    const fc2d_cudaclaw5_options_t* cudaclaw_options;
    cudaclaw_options = fc2d_cudaclaw5_get_options(glob);

    if (cudaclaw_options->ascii_out != 0)
    {
        fclaw2d_clawpatch_output_ascii(glob,iframe);
    }

    if (cudaclaw_options->vtk_out != 0)
    {
        fclaw2d_clawpatch_output_vtk(glob,iframe);
    }

}

/* ---------------------------------- Virtual table  ------------------------------------- */

static
fc2d_cudaclaw5_vtable_t* fc2d_cudaclaw5_vt_init()
{
    FCLAW_ASSERT(s_cudaclaw5_vt.is_set == 0);
    return &s_cudaclaw5_vt;
}


/* This is called from the user application. */
void fc2d_cudaclaw5_solver_initialize()
{
    int claw_version = 5;
    fclaw2d_clawpatch_vtable_initialize(claw_version);

    fclaw2d_vtable_t*          fclaw_vt = fclaw2d_vt();
    fclaw2d_patch_vtable_t*    patch_vt = fclaw2d_patch_vt();

    fc2d_cudaclaw5_vtable_t*   cuclaw5_vt = fc2d_cudaclaw5_vt_init();

    fclaw_vt->output_frame      = cudaclaw5_output;
    fclaw_vt->problem_setup     = cudaclaw5_setprob;    

    /* Default patch functions */
    patch_vt->initialize            = cudaclaw5_qinit;
    patch_vt->setup                 = cudaclaw5_setaux;
    patch_vt->physical_bc           = cudaclaw5_bc2;
    patch_vt->single_step_update    = cudaclaw5_update;

    cuclaw5_vt->b4step2   = cudaclaw5_b4step2;
    cuclaw5_vt->src2      = cudaclaw5_src2;

    /* Required functions  - error if NULL */
    cuclaw5_vt->fort_bc2       = CUDACLAW5_BC2_DEFAULT;
    cuclaw5_vt->fort_qinit     = NULL;
    cuclaw5_vt->fort_rpn2      = NULL;
    cuclaw5_vt->fort_rpt2      = NULL;

    /* Optional functions - call only if non-NULL */
    cuclaw5_vt->fort_setprob   = NULL;
    cuclaw5_vt->fort_setaux    = NULL;
    cuclaw5_vt->fort_b4step2   = NULL;
    cuclaw5_vt->fort_src2      = NULL;

    cuclaw5_vt->is_set = 1;
}


/* ----------------------------- User access to solver functions --------------------------- */


/* These are here in case the user wants to call Clawpack routines directly */

fc2d_cudaclaw5_vtable_t* fc2d_cudaclaw5_vt()
{
    FCLAW_ASSERT(s_cudaclaw5_vt.is_set != 0);
    return &s_cudaclaw5_vt;
}


void fc2d_cudaclaw5_setprob(fclaw2d_global_t *glob)
{
    cudaclaw5_setprob(glob);
}

void fc2d_cudaclaw5_setaux(fclaw2d_global_t *glob,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    cudaclaw5_setaux(glob,this_patch,this_block_idx,this_patch_idx);
}


/* This should only be called when a new fclaw2d_clawpatch_t is created. */
void fc2d_cudaclaw5_set_capacity(fclaw2d_global_t *glob,
                                  fclaw2d_patch_t *this_patch,
                                  int this_block_idx,
                                  int this_patch_idx)
{
    int mx,my,mbc,maux,mcapa;
    double dx,dy,xlower,ylower;
    double *aux, *area;
    fc2d_cudaclaw5_options_t *clawopt;

    clawopt = fc2d_cudaclaw5_get_options(glob);
    mcapa = clawopt->mcapa;

    fclaw2d_clawpatch_grid_data(glob,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    area = fclaw2d_clawpatch_get_area(glob,this_patch);

    fclaw2d_clawpatch_aux_data(glob,this_patch,&aux,&maux);
    FCLAW_ASSERT(maux >= mcapa && mcapa > 0);

    CUDACLAW5_SET_CAPACITY(&mx,&my,&mbc,&dx,&dy,area,&mcapa,
                            &maux,aux);
}

void fc2d_cudaclaw5_qinit(fclaw2d_global_t *glob,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    cudaclaw5_qinit(glob,this_patch,this_block_idx,this_patch_idx);
}

void fc2d_cudaclaw5_b4step2(fclaw2d_global_t* glob,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt)
{
    cudaclaw5_b4step2(glob,this_patch,this_block_idx,this_patch_idx,t,dt);
}

void fc2d_cudaclaw5_bc2(fclaw2d_global_t *glob,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt,
                         int intersects_bc[],
                         int time_interp)
{
    cudaclaw5_bc2(glob,this_patch,this_block_idx,this_block_idx,t,dt,
                   intersects_bc,time_interp);
}

void fc2d_cudaclaw5_src2(fclaw2d_global_t* glob,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    cudaclaw5_src2(glob,this_patch,this_block_idx,this_block_idx,t,dt);
}





