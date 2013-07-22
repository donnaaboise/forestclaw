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

#include "amr_options.h"

/* Proposed naming convention:
 * Parameter names in config file (= long option names) identical to the
 * C variable members of amr_options_t, except "-" in parameter name
 * corresponds to "_" in C variable.
 * For example the short option would be -F <filename> and the long option
 * --new-datafile=<Filename>.
 */

void
amr_options_add_int_array (sc_options_t * opt,
                           int opt_char, const char *opt_name,
                           const char **array_string,
                           const char *default_string,
                           int **int_array, int initial_length,
                           const char *help_string)
{
    *int_array = NULL;
    sc_options_add_string (opt, opt_char, opt_name,
                           array_string, default_string, help_string);
    amr_options_convert_int_array (*array_string, int_array, initial_length);
}

void
amr_options_convert_int_array (const char *array_string,
                               int **int_array, int new_length)
{
    int i;
    const char *beginptr;
    char *endptr;

    new_length = SC_MAX (new_length, 0);
    *int_array = SC_REALLOC (*int_array, int, new_length);

    beginptr = array_string;
    for (i = 0; i < new_length; ++i)
    {
        if (beginptr == NULL)
        {
            (*int_array)[i] = 0;
        }
        else
        {
            (*int_array)[i] = (int) strtol (beginptr, &endptr, 10);
            beginptr = endptr;
        }
    }
}

#if 0
/* This is now done in a post-processing step */
static void
amr_options_convert_arrays (amr_options_t * amropt)
{
    /*
    amr_options_convert_int_array (amropt->order_string, &amropt->order,
                                   SpaceDim);
    amr_options_convert_int_array (amropt->mthlim_string, &amropt->mthlim,
                                   amropt->mwaves);
     */
    amr_options_convert_int_array (amropt->mthbc_string, &amropt->mthbc,
                                   NumFaces);
}
#endif


amr_options_t *
amr_options_new (sc_options_t * opt)
{
    amr_options_t *amropt;

    amropt = SC_ALLOC_ZERO (amr_options_t, 1);

    sc_options_add_int (opt, 0, "mx", &amropt->mx, 8,
                        "Number of grid cells per patch in x");

    sc_options_add_int (opt, 0, "my", &amropt->my, 8,
                        "Number of grid cells per patch in y");

    sc_options_add_double (opt, 0, "initial_dt", &amropt->initial_dt, 0.1,
                           "Initial time step size");

    sc_options_add_int (opt, 0, "outstyle", &amropt->outstyle, 1,
                        "Output style (1,2,3)");

    /* If outstyle == 1 */
    sc_options_add_double (opt, 0, "tfinal", &amropt->tfinal, 1.0,
                           "Final time");

    sc_options_add_int (opt, 0, "nout", &amropt->nout, 10,
                        "Number of time steps");

    /* Only needed if outstyle == 3 */
    sc_options_add_int (opt, 0, "nstep", &amropt->nstep, 1,
                        "Number of steps to take between printing output files.");


    /* This is a hack to control the VTK output while still in development.
     * The values are numbers which can be bitwise-or'd together.
     * 0 - no VTK output ever.
     * 1 - output for all stages of amrinit.
     * 2 - output whenever amr_output() is called.
     */
    sc_options_add_int (opt, 0, "vtkout", &amropt->vtkout, 0,
                        "VTK output method");
    sc_options_add_double (opt, 0, "vtkspace", &amropt->vtkspace, 0.,
                           "VTK visual spacing");
    sc_options_add_int (opt, 0, "vtkwrite", &amropt->vtkwrite, 0,
                        "VTK write variant");

    /* output options */
    sc_options_add_int (opt, 0, "verbosity", &amropt->verbosity, 0,
                        "Verbosity mode [0]");
    sc_options_add_switch (opt, 0, "serialout", &amropt->serialout,
                           "Enable serial output [F]");
    sc_options_add_string (opt, 0, "prefix", &amropt->prefix, "fort",
                           "Output file prefix [fort]");

    /* more clawpack options */
    sc_options_add_double (opt, 0, "max_cfl", &amropt->max_cfl, 1,
                           "Maximum CFL allowed [1]");

    sc_options_add_double (opt, 0, "desired_cfl", &amropt->desired_cfl, 0.9,
                           "Maximum CFL allowed [0.9]");


    sc_options_add_int (opt, 0, "meqn", &amropt->meqn, 1,
                        "Number of equations [1]");

    sc_options_add_int (opt, 0, "mbc", &amropt->mbc, 2,
                        "Number of ghost cells [2]");

    /* Array of NumFaces many values */
    amr_options_add_int_array (opt, 0, "mthbc", &amropt->mthbc_string, NULL,
                               &amropt->mthbc, NumFaces,
                               "Physical boundary condition type");
    /* At this point amropt->mthbc is allocated. Set defaults if desired. */

    sc_options_add_int (opt, 0, "refratio", &amropt->refratio, 2,
                        "Refinement ratio (fixed) [2]");

    sc_options_add_int (opt, 0, "minlevel", &amropt->minlevel, 0,
                        "Minimum refinement level [0]");

    sc_options_add_int (opt, 0, "maxlevel", &amropt->maxlevel, 0,
                        "Maximum refinement level");

    sc_options_add_int (opt, 0, "regrid_interval", &amropt->regrid_interval,
                        0, "Regrid every ''regrid_interval'' steps");


    sc_options_add_double (opt, 0, "ax", &amropt->ax, 0, "xlower (ax)");
    sc_options_add_double (opt, 0, "bx", &amropt->bx, 1, "xupper (bx)");
    sc_options_add_double (opt, 0, "ay", &amropt->ay, 0, "ylower (ay)");
    sc_options_add_double (opt, 0, "by", &amropt->by, 1, "yupper (by)");

#if 0
    /* bool is not allocated, disable for now */
    /* Does bool get allocated somewhere? */
    sc_options_add_string (opt, 0, "manifold", &bool, "F", "Manifold [F]");

    amropt->manifold = bool[0] == 'T' ? 1 : 0;

    sc_options_add_string (opt, 0, "mapped", &bool, "F", "Mapped grid [F]");
    amropt->mapped = bool[0] == 'T' ? 1 : 0;
#endif

    /* ------------------------------------------------------------------- */
    /* Right now all switch options default to false, need to change that */
    /* I am okay with them being set to false by default, since I expect that
       user will set everything in the input file
     */
    sc_options_add_switch (opt, 0, "manifold", &amropt->manifold,
                           "Solution is on manifold [F]");
    sc_options_add_switch (opt, 0, "mapped", &amropt->mapped,
                           "Use mapped grid [F]");
    sc_options_add_switch (opt, 0, "use_fixed_dt", &amropt->use_fixed_dt,
                           "Use fixed coarse grid time step [F]");
    sc_options_add_switch (opt, 0, "check_conservation",
                           &amropt->check_conservation,
                           "Check conservation [F]");
    sc_options_add_switch (opt, 0, "subcycle", &amropt->subcycle,
                           "Use subcycling in time [F]");
    /* ------------------------------------------------------------------- */

    /* -----------------------------------------------------------------------
       Options will be read from this file, if a '-F' flag is used at the command
       line.  Use this file for local modifications that are not tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_add_inifile (opt, 'F', "inifile",
                            "Read waveprop options from this file");

    /* -----------------------------------------------------------------------
       This is the default file that will be read if no command line options are
       given.  This file is tracked by Git.
       ----------------------------------------------------------------------- */
    sc_options_load (sc_package_id, SC_LP_ALWAYS, opt,
                     "fclaw2d_defaults.ini");

    amr_postprocess_parms (amropt);

    return amropt;
}

void
amr_postprocess_parms (amr_options_t * amropt)
{
    /* -----------------------------------------------------------------------
       This has to happen after parameters have been parsed, but before we
       check parameters.
       ----------------------------------------------------------------------- */
    amr_options_convert_int_array (amropt->mthbc_string, &amropt->mthbc,
                                   NumFaces);

}


/* -----------------------------------------------------------------
   Check input parms
   ----------------------------------------------------------------- */
void
amr_checkparms (amr_options_t * gparms)
{
    /* Check outstyle. */
    if (gparms->outstyle == 1 && gparms->use_fixed_dt)
    {
        double dT_outer = gparms->tfinal / gparms->nout;
        double dT_inner = gparms->initial_dt;
        int nsteps = (int) floor (dT_outer / dT_inner + .5);
        if (fabs (nsteps * dT_inner - dT_outer) > 1e-8)
        {
            printf
                ("For fixed dt, initial time step size must divide tfinal/nout "
                 "exactly.\n");
            exit (1);
        }
    }

    /* Could also do basic sanity checks on mx,my,... */
}




#if 0
void
amr_options_parse (sc_options_t * opt, amr_options_t * amropt,
                   int argc, char **argv, int log_priority)
#endif
     void amr_options_parse (sc_options_t * opt, int argc, char **argv,
                             int log_priority)
{
    int retval;

    retval = sc_options_parse (sc_package_id, SC_LP_ERROR, opt, argc, argv);
    if (retval < 0)
    {
        sc_options_print_usage (sc_package_id, log_priority, opt, NULL);
        sc_abort_collective ("Option parsing failed");
    }
    sc_options_print_summary (sc_package_id, log_priority, opt);
#if 0
    if (sc_is_root ())
    {
        retval = sc_options_save (sc_package_id, SC_LP_ERROR, opt,
                                  "fclaw2d_defaults.ini.used");
        SC_CHECK_ABORT (!retval, "Option save failed");
    }
#endif

#if 0
    /* This is now done in a post-processing step */
    amr_options_convert_arrays (amropt);
#endif
}



void
amr_options_destroy (amr_options_t * amropt)
{
    /* These are now stored under amropt->waveprop_parms */
    /* SC_FREE (amropt->->order); */
    /* SC_FREE (amropt->mthlim); */
    SC_FREE (amropt->mthbc);
    SC_FREE (amropt);

}
