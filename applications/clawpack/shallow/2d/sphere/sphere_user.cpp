/*
Copyright (c) 2012-2023 Carsten Burstedde, Donna Calhoun
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

#include "sphere_user.h"

#include <fclaw_include_all.h>

#include <fclaw_clawpatch_pillow.h>
#include <fclaw_clawpatch.h>
#include <fclaw_clawpatch_options.h>

#include <fclaw2d_metric.h>

#include <fclaw2d_clawpatch_fort.h>

#include <fc2d_clawpack46.h> 
#include <fc2d_clawpack46_fort.h>
#include <fc2d_clawpack46_options.h>
#include <clawpack46_user_fort.h>


#include <fc2d_clawpack5.h> 
#include <fc2d_clawpack5_fort.h>
#include <fc2d_clawpack5_options.h>
#include <clawpack5_user_fort.h>

#if 0
// Fix syntax highlighting
#endif

/* This value needs to be   consistent with what is in setaux */
static int s_mbathy = 18;

static
void sphere_problem_setup(fclaw_global_t* glob)
{
    const user_options_t* user = sphere_get_options(glob);
    const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

    if (glob->mpirank == 0)
    {
        FILE *f = fopen("setprob_2d.data","w");
        fprintf(f,  "%-24d   %s",user->example,"\% example\n");
        fprintf(f,  "%-24.16f   %s",user->gravity,"\% gravity\n");
        fprintf(f,  "%-24d   %s",user->mapping,"\% mapping\n");
        fprintf(f,  "%-24d   %s",user->init_cond,"\% initial_condition\n");        

        fprintf(f,  "%-24.16f   %s",user->hin,"\% hin\n");
        fprintf(f,  "%-24.16f   %s",user->hout,"\% hout\n");

        fprintf(f,  "%-24.16f   %s",user->ring_inner,"\% ring-inner\n");
        fprintf(f,  "%-24.16f   %s",user->ring_outer,"\% ring-outer\n");
        fprintf(f,  "%-24d   %s",user->ring_units,"\% ring_units\n");

        fprintf(f,  "%-24.16f   %s",user->latitude[0],"\% latitude\n");
        fprintf(f,  "%-24.16f   %s",user->latitude[1],"\% latitude\n");
        fprintf(f,  "%-24.16f   %s",user->longitude[0],"\% longitude\n");
        fprintf(f,  "%-24.16f   %s",user->longitude[1],"\% longitude\n");

        fprintf(f,  "%-24.16f   %s",user->center[0],"\% center[0]\n");
        fprintf(f,  "%-24.16f   %s",user->center[1],"\% center[1]\n");
        fprintf(f,  "%-24.16f   %s",user->center[2],"\% center[2]\n");

        fprintf(f,  "%-24.16f   %s",user->theta_wave,"\% theta_wave\n");
        fprintf(f,  "%-24.16f   %s",user->theta_ridge,"\% theta_ridge\n");
        fprintf(f,  "%-24.16f   %s",user->gravity_ridge,"\% gravity_ridge\n");
        fprintf(f,  "%-24.16f   %s",user->ampl, "\% ampl\n");
        fprintf(f,  "%-24.16f   %s",user->alpha,"\% alpha\n");
        fprintf(f,  "%-24.16f   %s",user->speed,"\% speed\n");
        fprintf(f,  "%-24.16f   %s",user->bathy[0],"\% bathy[0]\n");
        fprintf(f,  "%-24.16f   %s",user->bathy[1],"\% bathy[1]\n");

        fprintf(f,  "%-24.16f   %s",fclaw_opt->refine_threshold,"\% refine_threshold\n");
        fprintf(f,  "%-24.16f   %s",fclaw_opt->coarsen_threshold,"\% coarsen_threshold\n");

        fprintf(f,  "%-24.16f   %s",fclaw_opt->phi,"\% phi\n");
        fprintf(f,  "%-24.16f   %s",fclaw_opt->theta,"\% theta\n");

        fclose(f);
    }
    fclaw_domain_barrier (glob->domain);
    SETPROB();
}

static
void sphere_patch_setup_manifold(fclaw_global_t *glob,
                                 fclaw_patch_t *patch,
                                 int blockno,
                                 int patchno)
{
    //const user_options_t *user = sphere_get_options(glob);

    int mx, my, mbc;
    double xlower, ylower, dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int maux;
    double *aux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    double *xd,*yd,*zd,*area;
    double *xp,*yp,*zp;
    fclaw_clawpatch_2d_metric_data(glob,patch,&xp,&yp,&zp,
                                  &xd,&yd,&zd,&area);

    double *xnormals, *ynormals,*xtangents,*ytangents;
    double *surfnormals, *edgelengths, *curvature;
    fclaw_clawpatch_2d_metric_data2(glob,patch,
                                   &xnormals,&ynormals,
                                   &xtangents,&ytangents,
                                   &surfnormals,
                                   &edgelengths,&curvature);


    // This is experimental ...
    //int* block_corner_count = fclaw_patch_block_corner_count(glob,patch);

    const user_options_t* user_opt = sphere_get_options(glob);
    if (user_opt->claw_version == 4)
    {
        CLAWPACK46_SET_BLOCK(&blockno);
        SPHERE_SETAUX(&mx,&my,&mbc,&xlower,&ylower,
                      &dx,&dy,area,xnormals,ynormals,
                      xtangents,ytangents,surfnormals, curvature,
                      edgelengths,
                      aux, &maux,&s_mbathy);
        CLAWPACK46_UNSET_BLOCK();
    }
    else
    {
        CLAWPACK5_SET_BLOCK(&blockno);
        SPHERE5_SETAUX(&mx,&my,&mbc,&xlower,&ylower,
                      &dx,&dy,area,xnormals,ynormals,
                      xtangents,ytangents,surfnormals, curvature,
                      edgelengths,
                      aux, &maux);
        CLAWPACK5_UNSET_BLOCK();
    }        
}

static
void cb_sphere_output_ascii(fclaw_domain_t *domain,
                             fclaw_patch_t *patch,
                             int blockno, int patchno,
                             void *user)
{
    fclaw_global_iterate_t* s = (fclaw_global_iterate_t*) user;
    fclaw_global_t *glob = (fclaw_global_t*) s->glob;

    int iframe = *((int *) s->user);    

    /* Get info not readily available to user */
    int local_num, global_num, level;
    fclaw_patch_get_info(glob->domain,patch,
                           blockno,patchno,
                           &global_num, 
                           &local_num,&level);

    int mx,my,mbc;
    double xlower,ylower,dx,dy;
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    double *q;
    int meqn;
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    double *aux;
    int maux;
    fclaw_clawpatch_aux_data(glob,patch,&aux,&maux);

    const user_options_t* user_opt = sphere_get_options(glob);
    if (user_opt->claw_version == 4)
        SPHERE_FORT_WRITE_FILE(&mx,&my,&meqn, &maux,&mbc,&xlower,&ylower,
                               &dx,&dy,q,aux,&iframe,&global_num,&level,
                               &blockno,&glob->mpirank);
    else
        SPHERE5_FORT_WRITE_FILE(&mx,&my,&meqn, &maux,&mbc,&xlower,&ylower,
                               &dx,&dy,q,aux,&iframe,&global_num,&level,
                               &blockno,&glob->mpirank);
}

static
void sphere_header_ascii(fclaw_global_t* glob,int iframe)
{
    double time = glob->curr_time;
    int ngrids = glob->domain->global_num_patches;

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int meqn = clawpatch_opt->meqn;
    int maux = clawpatch_opt->maux;

    SPHERE_FORT_WRITE_HEADER(&iframe,&time,&meqn,&maux,&ngrids);
}
/* ---------------------------------- AVERAGE ---------------------------------- */
static
void sphere_average_face(fclaw_global_t *glob,
                            fclaw_patch_t *coarse_patch,
                            fclaw_patch_t *fine_patch,
                            int idir,
                            int iface_coarse,
                            int p4est_refineFactor,
                            int refratio,
                            int time_interp,
                            int igrid,
                            fclaw_patch_transform_data_t* transform_data)
{
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);


    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    // const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);
    //fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);


    int maux;
    double *auxcoarse, *auxfine;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);    
    fclaw_clawpatch_aux_data(glob,fine_patch,   &auxfine,&maux);    

    int mx = clawpatch_opt->mx;
    int my = clawpatch_opt->my;

    fc2d_clawpack46_options_t *clawpack46_opt = 
                   fc2d_clawpack46_get_options(glob);
    int mcapa = clawpack46_opt->mcapa;

    /* These will be empty for non-manifolds cases */
    SPHERE_FORT_AVERAGE_FACE(&mx,&my,&mbc,&meqn,&mcapa, &s_mbathy, 
                             qcoarse,qfine,auxcoarse,auxfine, &maux, 
                             &idir,&iface_coarse, &igrid,
                             &transform_data);
}


static
void sphere_average_corner(fclaw_global_t *glob,
                              fclaw_patch_t *coarse_patch,
                              fclaw_patch_t *fine_patch,
                              int coarse_blockno,
                              int fine_blockno,
                              int coarse_corner,
                              int time_interp,
                              fclaw_patch_transform_data_t* transform_data)
{
    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    int maux;
    double *auxcoarse, *auxfine;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);    
    fclaw_clawpatch_aux_data(glob,fine_patch,   &auxfine,&maux);    

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    // this is experimental
    int* block_corner_count = fclaw_patch_block_corner_count(glob,coarse_patch);
    if (block_corner_count[coarse_corner] == 3)
        return;

    fc2d_clawpack46_options_t *clawpack46_opt = 
                   fc2d_clawpack46_get_options(glob);
    int mcapa = clawpack46_opt->mcapa;

    // const fclaw_options_t *fclaw_opt = fclaw_get_options(glob);
    // if (fill_ghost(glob,time_interp))
    {
        //fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;

        /* These will be empty for non-manifolds cases */
        SPHERE_FORT_AVERAGE_CORNER(&mx,&my,&mbc,&meqn,&mcapa,&s_mbathy,
                                   qcoarse,qfine,
                                   auxcoarse,auxfine,&maux,
                                   &coarse_corner,&transform_data);
    }
}

static
void sphere_average2coarse(fclaw_global_t *glob,
                           fclaw_patch_t *fine_patches,
                           fclaw_patch_t *coarse_patch,
                           int blockno, int fine0_patchno,
                           int coarse_patchno)

{
    const fclaw_clawpatch_options_t *clawpatch_opt = 
                    fclaw_clawpatch_get_options(glob);    
    int mbc = clawpatch_opt->mbc;
    int meqn = clawpatch_opt->meqn;


    for(int igrid = 0; igrid < fclaw_domain_num_siblings(glob->domain); igrid++)
    {
        fclaw_patch_t *fine_patch = &fine_patches[igrid];
        double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);
        double *qcoarse = fclaw_clawpatch_get_q(glob,coarse_patch);

        // const fclaw_options_t* fclaw_opt = fclaw_get_options(glob);

        // fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

        int maux;
        double *auxcoarse, *auxfine;
        fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);    
        fclaw_clawpatch_aux_data(glob,fine_patch,   &auxfine,&maux);    

        fc2d_clawpack46_options_t *clawpack46_opt = 
               fc2d_clawpack46_get_options(glob);
        int mcapa = clawpack46_opt->mcapa;


        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        SPHERE_FORT_AVERAGE2COARSE(&mx,&my,&mbc,&meqn,&mcapa, &s_mbathy, 
                                   qcoarse,qfine, auxcoarse, auxfine, 
                                   &maux, &igrid);
    }
}


/* ----------------------------------- INTERPOLATE ------------------------------------ */

static
void sphere_interpolate_face(fclaw_global_t *glob,
                             fclaw_patch_t *coarse_patch,
                             fclaw_patch_t *fine_patch,
                             int idir,
                             int iface_coarse,
                             int p4est_refineFactor,
                             int refratio,
                             int time_interp,
                             int igrid,
                             fclaw_patch_transform_data_t* transform_data)
{

    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);
    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    int mbc = clawpatch_opt->mbc;

    // Interpolation stencils will use corner data;  this allows us to 
    // set corner data at three-patch corners appropriately. 
    int* block_corner_count = fclaw_patch_block_corner_count(glob,coarse_patch);
    FCLAW2D_CLAWPATCH_SET_CORNER_COUNT(block_corner_count);

    int maux;
    double *auxcoarse, *auxfine;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);    
    fclaw_clawpatch_aux_data(glob,fine_patch,   &auxfine,&maux);    

    // if (fill_ghost(glob,time_interp))
    {
        //fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);

        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        SPHERE_FORT_INTERPOLATE_FACE(&mx,&my,&mbc,&meqn,
                                     qcoarse,qfine, auxcoarse, auxfine,
                                     &maux,&s_mbathy,&idir,
                                     &iface_coarse,&igrid,&transform_data);
    }
}

static
void sphere_interpolate_corner(fclaw_global_t* glob,
                                  fclaw_patch_t* coarse_patch,
                                  fclaw_patch_t* fine_patch,
                                  int coarse_blockno,
                                  int fine_blockno,
                                  int coarse_corner,
                                  int time_interp,
                                  fclaw_patch_transform_data_t* transform_data)

{
    const fclaw_clawpatch_options_t *clawpatch_opt = fclaw_clawpatch_get_options(glob);
    int mbc = clawpatch_opt->mbc;

    int meqn;
    double *qcoarse;
    fclaw_clawpatch_timesync_data(glob,coarse_patch,time_interp,&qcoarse,&meqn);

    double *qfine = fclaw_clawpatch_get_q(glob,fine_patch);

    // this is experimental
#if 0
    int* block_corner_count = fclaw_patch_block_corner_count(glob,coarse_patch);
    if (block_corner_count[coarse_corner] == 3)
        return;
#endif        


    int maux;
    double *auxcoarse, *auxfine;
    fclaw_clawpatch_aux_data(glob,coarse_patch,&auxcoarse,&maux);    
    fclaw_clawpatch_aux_data(glob,fine_patch,   &auxfine,&maux);    

    //if (fill_ghost(glob,time_interp))
    {
        //fclaw_clawpatch_vtable_t* clawpatch_vt = fclaw_clawpatch_vt(glob);
        int mx = clawpatch_opt->mx;
        int my = clawpatch_opt->my;
        SPHERE_FORT_INTERPOLATE_CORNER(&mx,&my,&mbc,&meqn,
                                       qcoarse,qfine,auxcoarse,auxfine,
                                       &maux,&s_mbathy,
                                       &coarse_corner,&transform_data);    
    }

}


static
void sphere_cons_check(fclaw_global_t *glob,
                       fclaw_patch_t *patch,
                       int blockno,
                       int patchno,
                       void *user)
{
    error_info_t* error_data = (error_info_t*) user;
    double* area = fclaw_clawpatch_get_2d_area(glob,patch);  /* Might be null */

    int meqn;
    double *q; 
    fclaw_clawpatch_soln_data(glob,patch,&q,&meqn);

    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);

    int mx, my, mbc;
    double xlower,ylower,dx,dy;

    FCLAW_ASSERT(clawpatch_vt->d2->fort_conservation_check != NULL);
    fclaw_clawpatch_2d_grid_data(glob,patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);
    SPHERE_FORT_CONSERVATION_CHECK(&blockno,&mx, &my, &mbc, 
                                   &meqn, &xlower, &ylower, &dx,&dy,
                                   area, q, error_data->mass,
                                   error_data->c_kahan);}


void sphere_link_solvers(fclaw_global_t *glob)
{
    /* ForestClaw core functions */
    fclaw_vtable_t *vt = fclaw_vt(glob);
    vt->problem_setup = &sphere_problem_setup;  /* Version-independent */

    fclaw_patch_vtable_t *patch_vt = fclaw_patch_vt(glob);
    patch_vt->setup   = &sphere_patch_setup_manifold;



    const user_options_t* user_opt = sphere_get_options(glob);
    if (user_opt->mapping == 3)
        fclaw_clawpatch_use_pillowsphere(glob);


    /* Clawpatch functions */    
    fclaw_clawpatch_vtable_t *clawpatch_vt = fclaw_clawpatch_vt(glob);
    clawpatch_vt->d2->fort_user_exceeds_threshold = &USER_EXCEEDS_THRESHOLD;

    clawpatch_vt->time_header_ascii = &sphere_header_ascii;
    clawpatch_vt->cb_output_ascii   = &cb_sphere_output_ascii;

    clawpatch_vt->conservation_check = sphere_cons_check;

#if 0
    /* This needs a C header */
    fclaw2d_metric_vtable_t *metric_vt = fclaw2d_metric_vt(glob);
    metric_vt->compute_area  = sphere_compute_area;
#endif    



    if (user_opt->claw_version == 4)
    {
        fc2d_clawpack46_vtable_t  *clawpack46_vt = fc2d_clawpack46_vt(glob);
        // clawpack46_vt->b4step2        = sphere_b4step2;
        clawpack46_vt->fort_qinit     = CLAWPACK46_QINIT;

        fc2d_clawpack46_options_t *clawpack46_opt = 
            fc2d_clawpack46_get_options(glob);

        if (clawpack46_opt->use_fwaves)
        {            
            clawpack46_vt->fort_rpn2 = &CLAWPACK46_RPN2_FWAVE; 
            clawpack46_vt->fort_rpt2 = &CLAWPACK46_RPT2_FWAVE;                  
        }
        else
        {
            clawpack46_vt->fort_rpn2 = &CLAWPACK46_RPN2;
            clawpack46_vt->fort_rpt2 = &CLAWPACK46_RPT2;            
        }
        clawpack46_vt->fort_rpn2_cons = &RPN2CONS_UPDATE_MANIFOLD;

        /* These are all 4.6 layout versions.  */
        patch_vt->average_face    = sphere_average_face;
        patch_vt->average_corner  = sphere_average_corner;
        patch_vt->average2coarse  = sphere_average2coarse;


        patch_vt->interpolate_face   = sphere_interpolate_face;
        patch_vt->interpolate_corner = sphere_interpolate_corner;
#if 0        
        clawpatch_vt->d2->fort_interpolate2fine   = SPHERE_FORT_INTERPOLATE2FINE;
#endif        

    }
    else
    {
        fc2d_clawpack5_vtable_t  *clawpack5_vt = fc2d_clawpack5_vt(glob);
        // clawpack46_vt->b4step2        = sphere_b4step2;
        clawpack5_vt->fort_qinit     = CLAWPACK5_QINIT;

        fc2d_clawpack5_options_t *clawpack5_opt = 
            fc2d_clawpack5_get_options(glob);
        if (clawpack5_opt->use_fwaves)
        {            
            clawpack5_vt->fort_rpn2 = &CLAWPACK5_RPN2_FWAVE; 
            clawpack5_vt->fort_rpt2 = &CLAWPACK5_RPT2_FWAVE;                  
        }
        else
        {
            clawpack5_vt->fort_rpn2 = &CLAWPACK5_RPN2;
            clawpack5_vt->fort_rpt2 = &CLAWPACK5_RPT2;            
        }

        clawpack5_vt->fort_rpn2_cons = &RPN2CONS_UPDATE_MANIFOLD;
    }
 }
