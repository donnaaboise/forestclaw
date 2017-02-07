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

#include <fclaw2d_global.h>

#include <fclaw2d_forestclaw.h>
#include <fclaw2d_clawpatch.hpp>

#include <fc2d_clawpack5.h>
#include "fc2d_geoclaw.h"
#include "fc2d_geoclaw_options.h"
#include "fclaw2d_map_brick.h"

/* Needed for debugging */
#include "types.h"

struct region_type region_type_for_debug;

static int s_geoclaw_package_id = -1;

static fc2d_geoclaw_vtable_t geoclaw_vt;

typedef struct fc2d_geoclaw_gauge_info
{
    sc_array_t *block_offsets;
    sc_array_t *coordinates;
} fc2d_geoclaw_gauge_info_t;

static fc2d_geoclaw_gauge_info_t gauge_info;

/* This is provided as a convencience to the user, and is
   called by the user app */
void fc2d_geoclaw_init_vtables(fclaw2d_vtable_t *fclaw_vt,
                               fc2d_geoclaw_vtable_t* geoclaw_vt)
{

    /* vt         : Functions required by ForestClaw
       geoclaw_vt : Specific to GeoClaw (or ClawPack) */

    fclaw2d_init_vtable(fclaw_vt);

    fclaw_vt->problem_setup      = &fc2d_geoclaw_setprob;  /* This function calls ... */
    geoclaw_vt->setprob          = NULL;                   /* ....     this function. */

    fclaw_vt->patch_setup        = &fc2d_geoclaw_patch_setup;
    geoclaw_vt->setaux           = &GEOCLAW_SETAUX;

    fclaw_vt->patch_initialize   = &fc2d_geoclaw_qinit;
    geoclaw_vt->qinit            = &GEOCLAW_QINIT;

    fclaw_vt->patch_physical_bc  = &fc2d_geoclaw_bc2;
    geoclaw_vt->bc2              = &GEOCLAW_BC2;

    fclaw_vt->patch_single_step_update = &fc2d_geoclaw_update;  /* Includes b4step2 and src2 */
    geoclaw_vt->b4step2          = &GEOCLAW_B4STEP2;
    geoclaw_vt->src2             = &GEOCLAW_SRC2;
    geoclaw_vt->rpn2             = &GEOCLAW_RPN2;
    geoclaw_vt->rpt2             = &GEOCLAW_RPT2;

    fclaw_vt->after_regrid             = &fc2d_geoclaw_after_regrid;
    fclaw_vt->regrid_tag4refinement    = &fc2d_geoclaw_patch_tag4refinement;
    fclaw_vt->regrid_tag4coarsening    = &fc2d_geoclaw_patch_tag4coarsening;

    fclaw_vt->regrid_interpolate2fine  = &fc2d_geoclaw_interpolate2fine;
    fclaw_vt->regrid_average2coarse    = &fc2d_geoclaw_average2coarse;

    fclaw_vt->write_header             = &fc2d_geoclaw_output_header_ascii;
    fclaw_vt->patch_write_file         = &fc2d_geoclaw_output_patch_ascii;

    /* diagnostic functions */
    fclaw_vt->fort_compute_error_norm = &FC2D_CLAWPACK5_FORT_COMPUTE_ERROR_NORM;
    fclaw_vt->fort_compute_patch_area = &FC2D_CLAWPACK5_FORT_COMPUTE_PATCH_AREA;
    fclaw_vt->fort_conservation_check = &FC2D_CLAWPACK5_FORT_CONSERVATION_CHECK;

    /* Patch functions */
    fclaw_vt->fort_copy_face     = &FC2D_CLAWPACK5_FORT_COPY_FACE;
    // fclaw_vt->fort_average_face = &FC2D_CLAWPACK5_FORT_AVERAGE_FACE;
    // fclaw_vt->fort_interpolate_face = &FC2D_CLAWPACK5_FORT_INTERPOLATE_FACE;

    fclaw_vt->average_face       = &fc2d_geoclaw_average_face;
    fclaw_vt->interpolate_face   = &fc2d_geoclaw_interpolate_face;

    fclaw_vt->fort_copy_corner   =  &FC2D_CLAWPACK5_FORT_COPY_CORNER;
    fclaw_vt->average_corner     =  &fc2d_geoclaw_average_corner;
    fclaw_vt->interpolate_corner =  &fc2d_geoclaw_interpolate_corner;
    // fclaw_vt->fort_average_corner     = &FC2D_CLAWPACK5_FORT_AVERAGE_CORNER;
    // fclaw_vt->fort_interpolate_corner = &FC2D_CLAWPACK5_FORT_INTERPOLATE_CORNER;

    /* Ghost Patch*/
    fclaw_vt->fort_ghostpack_qarea   = &FC2D_CLAWPACK5_FORT_GHOSTPACK_QAREA;
    fclaw_vt->ghostpatch_setup       = &fc2d_geoclaw_patch_setup;
    fclaw_vt->ghostpack_extra        = &fc2d_geoclaw_ghostpack_extra;

    fclaw_vt->fort_timeinterp = &FC2D_CLAWPACK5_FORT_TIMEINTERP;

    /* Update gauges */
    fclaw_vt->run_user_diagnostics = &fc2d_geoclaw_update_gauges;
}


void fc2d_geoclaw_set_vtables(fclaw2d_domain_t *domain,
                              fclaw2d_vtable_t *vt,
                              fc2d_geoclaw_vtable_t* geoclaw)
{
    geoclaw_vt = *geoclaw;
    fclaw2d_set_vtable(domain,vt);
}


/* Patch data is stored in each ClawPatch */
struct patch_aux_data
{
    FArrayBox auxarray;
    int maux;
};

fc2d_geoclaw_options_t* fc2d_geoclaw_get_options(fclaw2d_domain_t *domain)
{
    int id;
    fclaw_app_t* app;
    app = fclaw2d_domain_get_app(domain);
    id = fc2d_geoclaw_get_package_id();
    return (fc2d_geoclaw_options_t*) fclaw_package_get_options(app,id);
}

static patch_aux_data*
get_patch_data(ClawPatch *cp)
{
    patch_aux_data *wp =
        (patch_aux_data*) cp->clawpack_patch_data(s_geoclaw_package_id);
    return wp;
}

static void*
patch_data_new()
{
    patch_aux_data* data;
    data = new patch_aux_data;
    return (void*) data;
}

static void
patch_data_delete(void *data)
{
    patch_aux_data *pd = (patch_aux_data*) data;
    FCLAW_ASSERT(pd != NULL);
    delete pd;
}

static const fclaw_package_vtable_t geoclaw_patch_vtable = {
    patch_data_new,
    patch_data_delete
};


/* -----------------------------------------------------------
   Public interface to routines in this file
   ----------------------------------------------------------- */
void
fc2d_geoclaw_register_vtable (fclaw_package_container_t * pkg_container,
                                 fc2d_geoclaw_options_t * clawopt)
{
    FCLAW_ASSERT(s_geoclaw_package_id == -1);

    s_geoclaw_package_id =
      fclaw_package_container_add (pkg_container, clawopt,
                                   &geoclaw_patch_vtable);
}

void fc2d_geoclaw_register (fclaw_app_t* app, const char *configfile)
{
    fc2d_geoclaw_options_t* clawopt;
    int id;

    /* Register the options */
    clawopt = fc2d_geoclaw_options_register(app,configfile);

    /* And the package */

    /* Don't register a package more than once */
    FCLAW_ASSERT(s_geoclaw_package_id == -1);

    id = fclaw_package_container_add_pkg(app,clawopt,
                                         &geoclaw_patch_vtable);

    s_geoclaw_package_id = id;
}

int fc2d_geoclaw_get_package_id (void)
{
    return s_geoclaw_package_id;
}

static void
fc2d_geoclaw_define_auxarray_cp (fclaw2d_domain_t* domain, ClawPatch *cp)
{
    int maux;
    fc2d_geoclaw_options_t *clawpack_options;
    const amr_options_t* gparms;
    int mx,my,mbc;
    int ll[2], ur[2];

    gparms = get_domain_parms(domain);
    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    clawpack_options = fc2d_geoclaw_get_options(domain);
    maux = clawpack_options->maux;

    /* Construct an index box to hold the aux array */
    ll[0] = 1-mbc;
    ll[1] = 1-mbc;
    ur[0] = mx + mbc;
    ur[1] = my + mbc;    Box box(ll,ur);

    /* get solver specific data stored on this patch */
    patch_aux_data *clawpack_patch_data = get_patch_data(cp);
    clawpack_patch_data->auxarray.define(box,maux);
    clawpack_patch_data->maux = maux; // set maux in solver specific patch data (for completeness)
}

void fc2d_geoclaw_define_auxarray (fclaw2d_domain_t* domain,
                                      fclaw2d_patch_t* this_patch)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp (this_patch);
    fc2d_geoclaw_define_auxarray_cp (domain, cp);
}

static void
fc2d_geoclaw_get_auxarray_cp (fclaw2d_domain_t* domain,
                                 ClawPatch *cp, double **aux, int* maux)
{
    fc2d_geoclaw_options_t* clawpack_options;
    clawpack_options = fc2d_geoclaw_get_options(domain);
    *maux = clawpack_options->maux;

    patch_aux_data *clawpack_patch_data = get_patch_data(cp);
    *aux = clawpack_patch_data->auxarray.dataPtr();
}

void fc2d_geoclaw_aux_data(fclaw2d_domain_t* domain,
                              fclaw2d_patch_t *this_patch,
                              double **aux, int* maux)
{
    ClawPatch *cp = fclaw2d_clawpatch_get_cp (this_patch);
    fc2d_geoclaw_get_auxarray_cp (domain,cp,aux,maux);
}

void fc2d_geoclaw_maux(fclaw2d_domain_t* domain,int* maux)
{
    *maux = fc2d_geoclaw_get_maux(domain);
}

int fc2d_geoclaw_get_maux(fclaw2d_domain_t* domain)
{
    fc2d_geoclaw_options_t *clawpack_options;
    clawpack_options = fc2d_geoclaw_get_options(domain);
    return clawpack_options->maux;
}

void fc2d_geoclaw_setup(fclaw2d_domain_t *domain)
{
    const amr_options_t* gparms = get_domain_parms(domain);
    fc2d_geoclaw_options_t *geoclaw_options = fc2d_geoclaw_get_options(domain);

    GEOCLAW_SET_MODULES(&geoclaw_options->mwaves, &geoclaw_options->mcapa,
                        &gparms->meqn, &geoclaw_options->maux,
                        geoclaw_options->mthlim, geoclaw_options->method,
                        &gparms->ax, &gparms->bx, &gparms->ay, &gparms->by);
    fc2d_geoclaw_gauge_setup(domain);
}

void fc2d_geoclaw_gauge_setup(fclaw2d_domain_t* domain)
{
    const amr_options_t * gparms = get_domain_parms(domain);
    fc2d_geoclaw_options_t *geoclaw_options = fc2d_geoclaw_get_options(domain);

    /* --------------------------------------------------------
       Read gauges files 'gauges.data' to get number of gauges
       -------------------------------------------------------- */
    char fname[] = "gauges.data";
    int num = GEOCLAW_GAUGES_GETNUM(fname);
    int restart = 0;

    geoclaw_options->num_gauges = num;

    if (num == 0)
    {
        return;
    }

    geoclaw_options->gauges = FCLAW_ALLOC(geoclaw_gauge_t,num);

    /* Read gauges file for the locations, etc. of all gauges */
    GEOCLAW_GAUGES_INIT(&restart, &gparms->meqn, &num,  geoclaw_options->gauges, fname);

    /* -----------------------------------------------------
       Open gauge files and add header information
       ----------------------------------------------------- */
    char filename[14];    /* gaugexxxxx.txt */
    FILE *fp;

    for (int i = 0; i < geoclaw_options->num_gauges; ++i)
    {
        geoclaw_gauge_t g = geoclaw_options->gauges[i];
        sprintf(filename,"gauge%05d.txt",g.num);
        fp = fopen(filename, "w");
        fprintf(fp, "# gauge_id= %5d location=( %15.7e %15.7e ) num_eqn= %2d\n",
                g.num, g.xc, g.yc, gparms->meqn+1);
        fprintf(fp, "# Columns: level time h    hu    hv    eta\n");
        fclose(fp);
    }

    /* -----------------------------------------------------
       Set up block offsets and coordinate list for p4est
       search function
       ----------------------------------------------------- */

    fclaw2d_map_context_t* cont =
        fclaw2d_domain_get_map_context(domain);

    int is_brick = FCLAW2D_MAP_IS_BRICK(&cont);

    gauge_info.block_offsets = sc_array_new_size(sizeof(int), domain->num_blocks+1);
    gauge_info.coordinates = sc_array_new_size(2*sizeof(double), num);

    int *block_offsets = (int*) sc_array_index_int(gauge_info.block_offsets, 0);
    double *coordinates = (double*) sc_array_index_int(gauge_info.coordinates, 0);

    geoclaw_gauge_t *gauges = geoclaw_options->gauges;

    if (is_brick)
    {
        /* We don't know how the blocks are arranged in the brick domain
           so we reverse engineer this information
        */
        int nb,mi,mj,numblockgauge,numgaugeset;
        double x,y;
        double z;
        double xll,yll;
        double xur,yur;


        double x0,y0,x1,y1;
        x0 = 0;
        y0 = 0;
        x1 = 1;
        y1 = 1;

        FCLAW2D_MAP_BRICK_GET_DIM(&cont,&mi,&mj);

        numblockgauge = 0;
        numgaugeset = 0;
        block_offsets[0] = 0;
        for (nb = 0; nb < domain->num_blocks; ++nb)
        {
            /* Scale to [0,1]x[0,1], based on blockno */
            fclaw2d_map_c2m_nomap_brick(cont,nb,x0,y0,&xll,&yll,&z);
            fclaw2d_map_c2m_nomap_brick(cont,nb,x1,y1,&xur,&yur,&z);
            for(int i = 0; i < num; i++)
            {
                /* Map gauge to global [0,1]x[0,1] space */
                x = (gauges[i].xc - gparms->ax)/(gparms->bx-gparms->ax);
                y = (gauges[i].yc - gparms->ay)/(gparms->by-gparms->ay);
                if (xll <= x && x <= xur && yll <= y && y <= yur)
                {
                    gauges[i].blockno = nb;
                    coordinates[2*numgaugeset] = mi*(x - xll);
                    coordinates[2*numgaugeset+1] = mj*(y - yll);
                    gauges[i].location_in_results = numgaugeset;
                    numblockgauge++;
                    numgaugeset++;
                }
            }
            block_offsets[nb+1] = numblockgauge;
            numblockgauge = 0;
        }
    }
    else
    {
        for(int i = 0; i < num; i++)
        {
            gauges[i].blockno = 0;
            gauges[i].location_in_results = i;
        }

        block_offsets[0] = 0;
        block_offsets[1] = num;

        for (int i = 0; i < num; ++i)
        {
            coordinates[2*i] = (gauges[i].xc - gparms->ax)/(gparms->bx-gparms->ax);
            coordinates[2*i+1] = (gauges[i].yc - gparms->ay)/(gparms->by-gparms->ay);
        }
    }
}


void fc2d_geoclaw_update_gauges(fclaw2d_domain_t *domain, const double tcurr)
{
    int mx,my,mbc,meqn,maux;
    double dx,dy,xlower,ylower,eta;
    double *q, *aux;
    double *var;
    char filename[14];  /* gaugexxxxx.txt */
    FILE *fp;

    const fc2d_geoclaw_options_t *geoclaw_options;
    const amr_options_t* gparms;
    geoclaw_options = fc2d_geoclaw_get_options(domain);
    gparms = get_domain_parms(domain);

    if (geoclaw_options->num_gauges == 0)
    {
        return;
    }

    fclaw2d_block_t *block;
    fclaw2d_patch_t *patch;

    var = FCLAW_ALLOC(double, gparms->meqn);
    geoclaw_gauge_t *gauges = geoclaw_options->gauges;
    geoclaw_gauge_t g;
    for (int i = 0; i < geoclaw_options->num_gauges; ++i)
    {
        g = gauges[i];
        block = &domain->blocks[g.blockno];
        if (g.patchno >= 0)
        {
            patch = &block->patches[g.patchno];
            FCLAW_ASSERT(patch != NULL);
            fclaw2d_clawpatch_grid_data(domain,patch,&mx,&my,&mbc,
                                        &xlower,&ylower,&dx,&dy);

            fclaw2d_clawpatch_soln_data(domain,patch,&q,&meqn);
            fc2d_geoclaw_aux_data(domain,patch,&aux,&maux);

            FCLAW_ASSERT(g.xc >= xlower && g.xc <= xlower+mx*dx);
            FCLAW_ASSERT(g.yc >= ylower && g.yc <= ylower+my*dy);
            if (tcurr >= g.t1 && tcurr <= g.t2)
            {
                GEOCLAW_UPDATE_GAUGE(&mx,&my,&mbc,&meqn,&xlower,&ylower,&dx,&dy,
                                      q,&maux,aux,&g.xc,&g.yc,var,&eta);
                sprintf(filename,"gauge%05d.txt",g.num);
                fp = fopen(filename, "a");
                fprintf(fp, "%5d %15.7e %15.7e %15.7e %15.7e %15.7e\n",
                        patch->level,tcurr,var[0],var[1],var[2],eta);
                fclose(fp);
            }
        }
    }
    FCLAW_FREE(var);
}

void fc2d_geoclaw_finalize(fclaw2d_domain_t *domain)
{
    sc_array_destroy(gauge_info.block_offsets);
    sc_array_destroy(gauge_info.coordinates);
}


void fc2d_geoclaw_after_regrid(fclaw2d_domain_t *domain)
{
    int i,index;

    fc2d_geoclaw_options_t *geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    /* Locate each gauge in the new mesh */
    int num = geoclaw_options->num_gauges;

    if (num == 0)
    {
        return;
    }

    sc_array_t *results = sc_array_new_size(sizeof(int), num);
    fclaw2d_domain_search_points(domain, gauge_info.block_offsets,
                                 gauge_info.coordinates, results);

    for (i = 0; i < geoclaw_options->num_gauges; ++i)
    {
        index = geoclaw_options->gauges[i].location_in_results;
        FCLAW_ASSERT(index >= 0 && index < num);

        /* patchno == -1  : Patch is not on this processor
           patchno >= 0   : Patch number is local patch list.
        */

        int patchno = *((int *) sc_array_index_int(results, index));
        geoclaw_options->gauges[i].patchno = patchno;
    }
    sc_array_destroy(results);
}


void fc2d_geoclaw_setprob(fclaw2d_domain_t *domain)
{
    if (geoclaw_vt.setprob != NULL)
    {
        geoclaw_vt.setprob();
    }
}


void fc2d_geoclaw_patch_setup(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx)
{
    fc2d_geoclaw_setaux(domain,this_patch,this_block_idx,this_patch_idx);
}

void fc2d_geoclaw_ghostpatch_setup(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *this_patch,
                                   int this_block_idx,
                                   int this_patch_idx)
{
    fc2d_geoclaw_options_t* geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    if(geoclaw_options->ghost_patch_pack_aux == 0){
        fc2d_geoclaw_setaux(domain,this_patch,this_block_idx,this_patch_idx);
    }
}

/* This should only be called when a new ClawPatch is created. */
void fc2d_geoclaw_setaux(fclaw2d_domain_t *domain,
                            fclaw2d_patch_t *this_patch,
                            int this_block_idx,
                            int this_patch_idx)
{
    FCLAW_ASSERT(geoclaw_vt.setaux != NULL);
    int mx,my,mbc,maux;
    double xlower,ylower,dx,dy;
    double *aux;

    int is_ghost = fclaw2d_patch_is_ghost(this_patch);

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    int mint = 2*mbc;
    int nghost = mbc;

    fc2d_geoclaw_define_auxarray (domain,this_patch);

    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);
    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt.setaux(&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,
                      &maux,aux,&is_ghost,&nghost,&mint);
    GEOCLAW_UNSET_BLOCK();
}


void fc2d_geoclaw_qinit(fclaw2d_domain_t *domain,
                           fclaw2d_patch_t *this_patch,
                           int this_block_idx,
                           int this_patch_idx)
{
    FCLAW_ASSERT(geoclaw_vt.qinit != NULL); /* Must initialized */
    int mx,my,mbc,meqn,maux;
    double dx,dy,xlower,ylower;
    double *q, *aux;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    /* Call to classic Clawpack 'qinit' routine.  This must be user defined */
    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt.qinit(&meqn,&mbc,&mx,&my,&xlower,&ylower,&dx,&dy,q,
                     &maux,aux);
    GEOCLAW_UNSET_BLOCK();
}

void fc2d_geoclaw_b4step2(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t, double dt)

{
    FCLAW_ASSERT(geoclaw_vt.b4step2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt.b4step2(&mbc,&mx,&my,&meqn,q,&xlower,&ylower,
                       &dx,&dy,&t,&dt,&maux,aux);
    GEOCLAW_UNSET_BLOCK();
}

void fc2d_geoclaw_src2(fclaw2d_domain_t *domain,
                          fclaw2d_patch_t *this_patch,
                          int this_block_idx,
                          int this_patch_idx,
                          double t,
                          double dt)
{
    FCLAW_ASSERT(geoclaw_vt.src2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    GEOCLAW_SET_BLOCK(&this_block_idx);
    geoclaw_vt.src2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                    &dx,&dy,q,&maux,aux,&t,&dt);
    GEOCLAW_UNSET_BLOCK();

}



/* Use this to return only the right hand side of the clawpack algorithm */
double fc2d_geoclaw_step2_rhs(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *this_patch,
                                 int this_block_idx,
                                 int this_patch_idx,
                                 double t,
                                 double *rhs)
{
    /* This should evaluate the right hand side, but not actually do the update.
       This will be useful in cases where we want to use something other than
       a single step method.  For example, in a RK scheme, one might want to
       call the right hand side to evaluate stages. */
    return 0;
}


void fc2d_geoclaw_bc2(fclaw2d_domain *domain,
                         fclaw2d_patch_t *this_patch,
                         int this_block_idx,
                         int this_patch_idx,
                         double t,
                         double dt,
                         fclaw_bool intersects_phys_bdry[],
                         fclaw_bool time_interp)
{
    FCLAW_ASSERT(geoclaw_vt.bc2 != NULL);

    int mx,my,mbc,meqn, maux;
    double xlower,ylower,dx,dy;
    double *aux,*q;

    fclaw2d_clawpatch_grid_data(domain,this_patch, &mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    fclaw2d_block_t *this_block = &domain->blocks[this_block_idx];
    fclaw2d_block_data_t *bdata = fclaw2d_block_get_data(this_block);
    int *block_mthbc = bdata->mthbc;

    /* Set a local copy of mthbc that can be used for a patch. */
    int mthbc[NumFaces];
    for(int i = 0; i < NumFaces; i++)
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
    fclaw2d_clawpatch_timesync_data(domain,this_patch,time_interp,&q,&meqn);

    GEOCLAW_SET_BLOCK(&this_block_idx);
    /*
    geoclaw_vt.bc2(q, aux, nrow, ncol ,meqn, maux, hx, hy,
                   level, time, xleft, xright, ybot, ytop,
                   &xlower, &ylower, xupper, yupper, 0,0,0);
    */
    geoclaw_vt.bc2(&meqn,&mbc,&mx,&my,&xlower,&ylower,
                   &dx,&dy,q,&maux,aux,&t,&dt,mthbc);
    GEOCLAW_UNSET_BLOCK();

}


/* This is called from the single_step callback. and is of type 'flaw_single_step_t' */
double fc2d_geoclaw_step2(fclaw2d_domain_t *domain,
                             fclaw2d_patch_t *this_patch,
                             int this_block_idx,
                             int this_patch_idx,
                             double t,
                             double dt)
{
    fc2d_geoclaw_options_t* geoclaw_options;
    ClawPatch *cp;
    int level;
    double *qold, *aux;
    int mx, my, meqn, maux, mbc;
    double xlower, ylower, dx,dy;

    FCLAW_ASSERT(geoclaw_vt.rpn2 != NULL);
    FCLAW_ASSERT(geoclaw_vt.rpt2 != NULL);

    geoclaw_options = fc2d_geoclaw_get_options(domain);

    cp = fclaw2d_clawpatch_get_cp(this_patch);

    level = this_patch->level;

    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    cp->save_current_step();  // Save for time interpolation

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&qold,&meqn);

    int mwaves = geoclaw_options->mwaves;

    int maxm = fmax(mx,my);

    double cflgrid = 0.0;

    /* This is no longer needed */
    int mwork = (maxm+2*mbc)*(12*meqn + (meqn+1)*mwaves + 3*maux + 2);
    double* work = new double[mwork];

    int size = meqn*(mx+2*mbc)*(my+2*mbc);
    double* fp = new double[size];
    double* fm = new double[size];
    double* gp = new double[size];
    double* gm = new double[size];

    int* block_corner_count = fclaw2d_patch_block_corner_count(domain,this_patch);

    GEOCLAW_STEP2_WRAP(&maxm, &meqn, &maux, &mbc, geoclaw_options->method,
                       geoclaw_options->mthlim, &geoclaw_options->mcapa,
                       &mwaves,&mx, &my, qold, aux, &dx, &dy, &dt, &cflgrid,
                       work, &mwork, &xlower, &ylower, &level,&t, fp, fm, gp, gm,
                       geoclaw_vt.rpn2, geoclaw_vt.rpt2,
                       block_corner_count);

    /* Accumulate fluxes needed for conservative fix-up */
    if (geoclaw_vt.fluxfun != NULL)
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

double fc2d_geoclaw_update(fclaw2d_domain_t *domain,
                              fclaw2d_patch_t *this_patch,
                              int this_block_idx,
                              int this_patch_idx,
                              double t,
                              double dt)
{
    const fc2d_geoclaw_options_t* geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    GEOCLAW_TOPO_UPDATE(&t);

    if (geoclaw_vt.b4step2 != NULL)
    {
        fc2d_geoclaw_b4step2(domain,
                             this_patch,
                             this_block_idx,
                             this_patch_idx,t,dt);
    }

    double maxcfl = fc2d_geoclaw_step2(domain,
                                       this_patch,
                                       this_block_idx,
                                       this_patch_idx,t,dt);
    if (geoclaw_options->src_term > 0)
    {
        fc2d_geoclaw_src2(domain,
                          this_patch,
                          this_block_idx,
                          this_patch_idx,t,dt);
    }

    return maxcfl;
}

int fc2d_geoclaw_patch_tag4refinement(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *this_patch,
                                      int blockno, int this_patch_idx,
                                      int initflag)
{
    int mx,my,mbc,meqn,maux,mbathy;
    double xlower,ylower,dx,dy, t;
    double *q, *aux;
    int tag_patch;
    int level,maxlevel;
    const amr_options_t * gparms = get_domain_parms(domain);
    const fc2d_geoclaw_options_t* geoclaw_options;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);
    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    level = this_patch->level;
    maxlevel = gparms->maxlevel;
    t = fclaw2d_domain_get_time(domain);

    geoclaw_options = fc2d_geoclaw_get_options(domain);
    mbathy = geoclaw_options->mbathy;

    tag_patch = 0;

    FC2D_GEOCLAW_FORT_TAG4REFINEMENT(&mx,&my,&mbc,&meqn,&maux,&xlower,&ylower,
                                     &dx,&dy,&t,&blockno,q,aux,&mbathy,&level,&maxlevel,
                                     &initflag,&tag_patch);

    return tag_patch;
}


int fc2d_geoclaw_patch_tag4coarsening(fclaw2d_domain_t *domain,
                                      fclaw2d_patch_t *fine_patches,
                                      int blockno, int patchno)

{
    fclaw2d_vtable_t vt;
    const amr_options_t *amropt;
    fc2d_geoclaw_options_t* geoclaw_options;

    int mx,my,mbc,meqn,maux,mbathy;
    double xlower,ylower,dx,dy,t;
    double *q[4], *aux[4];
    int tag_patch,igrid;
    //double coarsen_threshold;
    int level, maxlevel;

    fclaw2d_patch_t *patch0;
    patch0 = &fine_patches[0];


    amropt = get_domain_parms(domain);
    geoclaw_options = fc2d_geoclaw_get_options(domain);
    mbathy = geoclaw_options->mbathy;
    //coarsen_threshold = amropt->coarsen_threshold;
    t = fclaw2d_domain_get_time(domain);

    vt = fclaw2d_get_vtable(domain);

    fclaw2d_clawpatch_grid_data(domain,patch0,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    for (igrid = 0; igrid < 4; igrid++)
    {
        fclaw2d_clawpatch_soln_data(domain,&fine_patches[igrid],&q[igrid],&meqn);
        fc2d_geoclaw_aux_data(domain,&fine_patches[igrid],&aux[igrid],&maux);
    }

    level = fine_patches[0].level;
    maxlevel = amropt->maxlevel;

    tag_patch = 0;
    FC2D_GEOCLAW_FORT_TAG4COARSENING(&patchno,&mx,&my,&mbc,&meqn,&maux,&xlower,&ylower,&dx,&dy,
                                     &t,q[0],q[1],q[2],q[3],aux[0],aux[1],aux[2],aux[3],
                                     &mbathy,&level,&maxlevel, &geoclaw_options->dry_tolerance_c,
                                     &geoclaw_options->wave_tolerance_c,
                                     &geoclaw_options->speed_tolerance_entries_c,
                                     geoclaw_options->speed_tolerance_c, &tag_patch);
    /* Print patches level, tag_patch, xlower, ylower to observe the refinement*/
    // fclaw_global_infof("tag_patch %d, xlower %f, ylower %f, level %d\n",
    //                     tag_patch, xlower, ylower, level);
    return tag_patch;

}

void fc2d_geoclaw_interpolate2fine(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *coarse_patch,
                                   fclaw2d_patch_t *fine_patches,
                                   int this_blockno, int coarse_patchno,
                                   int fine0_patchno)

{
    fclaw2d_vtable_t vt;
    int mx,my,mbc,meqn,maux,refratio,p4est_refineFactor,mbathy;
    double *qcoarse,*qfine,*auxcoarse,*auxfine;
    // double *areacoarse,*areafine;
    // double *xp,*yp,*zp,*xd,*yd,*zd;
    int igrid;

    const amr_options_t* gparms;
    const fc2d_geoclaw_options_t*  geoclaw_options;
    fclaw2d_patch_t* fine_patch;

    vt = fclaw2d_get_vtable(domain);

    gparms = get_domain_parms(domain);
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    mbathy = geoclaw_options->mbathy;
    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;
    refratio = gparms->refratio;
    p4est_refineFactor = FCLAW2D_P4EST_REFINE_FACTOR;

    // fclaw2d_clawpatch_metric_data(domain,coarse_patch,&xp,&yp,&zp,
    //                               &xd,&yd,&zd,&areacoarse);

    fclaw2d_clawpatch_soln_data(domain,coarse_patch,&qcoarse,&meqn);
    fc2d_geoclaw_aux_data(domain,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw2d_clawpatch_soln_data(domain,fine_patch,&qfine,&meqn);
        fc2d_geoclaw_aux_data(domain,fine_patch,&auxfine,&maux);

        // if (gparms->manifold)
        // {
        //     fclaw2d_clawpatch_metric_data(domain,fine_patch,&xp,&yp,&zp,
        //                                   &xd,&yd,&zd,&areafine);
        // }
        FC2D_GEOCLAW_FORT_INTERPOLATE2FINE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                           &maux,auxcoarse,auxfine,&mbathy,
                                           &p4est_refineFactor,&refratio,
                                           &igrid);
        // vt.fort_interpolate2fine(&mx,&my,&mbc,&meqn,qcoarse,qfine,
        //                          areacoarse, areafine, &igrid,
        //                          &gparms->manifold);

    }
}

void fc2d_geoclaw_average2coarse(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *fine_patches,
                                 fclaw2d_patch_t *coarse_patch,
                                 int blockno, int fine0_patchno,
                                 int coarse_patchno)

{
    fclaw2d_vtable_t vt;
    int mx,my,mbc,meqn,maux,refratio,p4est_refineFactor,mbathy;
    double *qcoarse,*qfine,*auxcoarse,*auxfine;
    // double *areacoarse,*areafine;
    // double *xp,*yp,*zp,*xd,*yd,*zd;
    int igrid,mcapa;

    const amr_options_t* gparms;
    const fc2d_geoclaw_options_t*  geoclaw_options;

    fclaw2d_patch_t* fine_patch;

    vt = fclaw2d_get_vtable(domain);

    gparms = get_domain_parms(domain);
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    mcapa = geoclaw_options->mcapa;
    mx  = gparms->mx;
    my  = gparms->my;
    mbc = gparms->mbc;
    refratio = gparms->refratio;
    p4est_refineFactor = FCLAW2D_P4EST_REFINE_FACTOR;
    mbathy = geoclaw_options->mbathy;

    // fclaw2d_clawpatch_metric_data(domain,coarse_patch,&xp,&yp,&zp,
    //                               &xd,&yd,&zd,&areacoarse);

    fclaw2d_clawpatch_soln_data(domain,coarse_patch,&qcoarse,&meqn);
    fc2d_geoclaw_aux_data(domain,coarse_patch,&auxcoarse,&maux);

    /* Loop over four siblings (z-ordering) */
    for (igrid = 0; igrid < 4; igrid++)
    {
        fine_patch = &fine_patches[igrid];

        fclaw2d_clawpatch_soln_data(domain,fine_patch,&qfine,&meqn);
        fc2d_geoclaw_aux_data(domain,fine_patch,&auxfine,&maux);

        // if (gparms->manifold)
        // {
        //     fclaw2d_clawpatch_metric_data(domain,fine_patch,&xp,&yp,&zp,
        //                                   &xd,&yd,&zd,&areafine);
        // }
        FC2D_GEOCLAW_FORT_AVERAGE2COARSE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                         &maux,auxcoarse,auxfine,&mcapa,&mbathy,
                                         &p4est_refineFactor,&refratio,
                                         &igrid);

    }
}

void fc2d_geoclaw_average_face(fclaw2d_domain_t *domain,
                                    fclaw2d_patch_t *coarse_patch,
                                    fclaw2d_patch_t *fine_patch,
                                    int idir,
                                    int iface_coarse,
                                    int p4est_refineFactor,
                                    int refratio,
                                    fclaw_bool time_interp,
                                    int igrid,
                                    fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mbc;
    int mcapa, mbathy, maux;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;
    const amr_options_t *gparms = get_domain_parms(domain);
    const fc2d_geoclaw_options_t *geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    /* These will be empty for non-manifolds cases */
    fc2d_geoclaw_aux_data(domain,coarse_patch,&auxcoarse,&maux);
    fc2d_geoclaw_aux_data(domain,fine_patch,&auxfine,&maux);

    mcapa = geoclaw_options->mcapa;
    mbathy = geoclaw_options->mbathy;

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    int manifold = gparms->manifold;
    FC2D_GEOCLAW_FORT_AVERAGE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,
                                   &maux,auxcoarse,auxfine,&mcapa,&mbathy,
                                   &idir,&iface_coarse,&p4est_refineFactor,&refratio,
                                   &igrid,&manifold,&transform_data);
}

void fc2d_geoclaw_interpolate_face(fclaw2d_domain_t *domain,
                                   fclaw2d_patch_t *coarse_patch,
                                   fclaw2d_patch_t *fine_patch,
                                   int idir,
                                   int iside,
                                   int p4est_refineFactor,
                                   int refratio,
                                   fclaw_bool time_interp,
                                   int igrid,
                                   fclaw2d_transform_data_t* transform_data)
{

    int meqn,mx,my,mbc,maux,mbathy;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;
    const amr_options_t *gparms = get_domain_parms(domain);
    const fc2d_geoclaw_options_t *geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    fc2d_geoclaw_aux_data(domain,coarse_patch,&auxcoarse,&maux);
    fc2d_geoclaw_aux_data(domain,fine_patch,&auxfine,&maux);

    mbathy = geoclaw_options->mbathy;
    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    FC2D_GEOCLAW_FORT_INTERPOLATE_FACE(&mx,&my,&mbc,&meqn,qcoarse,qfine,&maux,
                                       auxcoarse,auxfine,&mbathy,
                                       &idir, &iside,
                                       &p4est_refineFactor,&refratio,&igrid,
                                       &transform_data);
}

void fc2d_geoclaw_average_corner(fclaw2d_domain_t *domain,
                                 fclaw2d_patch_t *coarse_patch,
                                 fclaw2d_patch_t *fine_patch,
                                 int coarse_corner,
                                 int refratio,
                                 fclaw_bool time_interp,
                                 fclaw2d_transform_data_t* transform_data)
{
    int meqn,mx,my,mbc;
    int maux,mbathy,mcapa;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;
    const amr_options_t *gparms = get_domain_parms(domain);
    const fc2d_geoclaw_options_t *geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    fc2d_geoclaw_aux_data(domain,coarse_patch,&auxcoarse,&maux);
    fc2d_geoclaw_aux_data(domain,fine_patch,&auxfine,&maux);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;

    mcapa=geoclaw_options->mcapa;
    mbathy=geoclaw_options->mbathy;

    int manifold = gparms->manifold;
    FC2D_GEOCLAW_FORT_AVERAGE_CORNER(&mx,&my,&mbc,&meqn,&refratio,
                                     qcoarse,qfine,&maux,auxcoarse,auxfine,
                                     &mcapa,&mbathy,&manifold,&coarse_corner,
                                     &transform_data);

}

void fc2d_geoclaw_interpolate_corner(fclaw2d_domain_t* domain,
                                     fclaw2d_patch_t* coarse_patch,
                                     fclaw2d_patch_t* fine_patch,
                                     int coarse_corner,
                                     int refratio,
                                     fclaw_bool time_interp,
                                     fclaw2d_transform_data_t* transform_data)

{
    int meqn,mx,my,mbc;
    int mbathy,maux;
    double *qcoarse, *qfine;
    double *auxcoarse, *auxfine;

    const amr_options_t *gparms = get_domain_parms(domain);
    const fc2d_geoclaw_options_t *geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    mx = gparms->mx;
    my = gparms->my;
    mbc = gparms->mbc;
    meqn = gparms->meqn;

    mbathy = geoclaw_options->mbathy;

    fclaw2d_clawpatch_timesync_data(domain,coarse_patch,time_interp,&qcoarse,&meqn);
    qfine = fclaw2d_clawpatch_get_q(domain,fine_patch);

    fc2d_geoclaw_aux_data(domain,coarse_patch,&auxcoarse,&maux);
    fc2d_geoclaw_aux_data(domain,fine_patch,&auxfine,&maux);

    FC2D_GEOCLAW_FORT_INTERPOLATE_CORNER(&mx,&my,&mbc,&meqn,
                                         &refratio,qcoarse,qfine,&maux,
                                         auxcoarse,auxfine,&mbathy,
                                         &coarse_corner,&transform_data);

}

void fc2d_geoclaw_output_header_ascii(fclaw2d_domain_t* domain,
                                      int iframe)
{
    const amr_options_t *amropt;
    const fc2d_geoclaw_options_t *geoclaw_options;
    int meqn,maux,ngrids;
    double time;

    amropt = fclaw2d_forestclaw_get_options(domain);
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    time = fclaw2d_domain_get_time(domain);
    ngrids = fclaw2d_domain_get_num_patches(domain);

    meqn = amropt->meqn;
    maux = geoclaw_options->maux;

    FC2D_GEOCLAW_FORT_WRITE_HEADER(&iframe,&time,&meqn,&maux,&ngrids);

    /* Is this really necessary? */
    /* FCLAW2D_OUTPUT_NEW_QFILE(&iframe); */
}

void fc2d_geoclaw_output_patch_ascii(fclaw2d_domain_t *domain,
                                     fclaw2d_patch_t *this_patch,
                                     int this_block_idx, int this_patch_idx,
                                     int iframe,int patch_num,int level)
{
    const fc2d_geoclaw_options_t *geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

    int mx,my,mbc,meqn,maux;
    int mbathy=geoclaw_options->mbathy;
    double xlower,ylower,dx,dy;
    double *q,*aux;

    fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                &xlower,&ylower,&dx,&dy);

    fclaw2d_clawpatch_soln_data(domain,this_patch,&q,&meqn);
    fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);

    FC2D_GEOCLAW_FORT_WRITE_FILE(&mx,&my,&meqn,&maux,&mbathy,&mbc,&xlower,&ylower,
                                 &dx,&dy,q,aux,&iframe,&patch_num,&level,
                                 &this_block_idx,&domain->mpirank);
}

void fc2d_geoclaw_ghostpack_extra(fclaw2d_domain_t *domain,
                                  fclaw2d_patch_t *this_patch,
                                  int mint,
                                  double *auxpack,
                                  int auxsize, int packmode, 
                                  int* ierror)
{
    const fc2d_geoclaw_options_t *geoclaw_options;
    geoclaw_options = fc2d_geoclaw_get_options(domain);

  
    if (geoclaw_options->ghost_patch_pack_aux)
    {
        int mx,my,mbc,maux;
        double xlower,ylower,dx,dy;
        double *aux;

        fclaw2d_clawpatch_grid_data(domain,this_patch,&mx,&my,&mbc,
                                    &xlower,&ylower,&dx,&dy);
        fc2d_geoclaw_aux_data(domain,this_patch,&aux,&maux);
        FC2D_GEOCLAW_FORT_GHOSTPACKAUX(&mx,&my,&mbc,&maux,
                                   &mint,aux,auxpack,&auxsize,
                                   &packmode,ierror);  
    }
}