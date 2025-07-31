subroutine setprob()
    implicit none

    double precision pi, pi2, deg2rad
    common /compi/ pi, pi2, deg2rad

    integer example
    common /common_ex/ example

    double precision revs_per_second
    common /com_latlong/ revs_per_second

    integer init_cond
    common /latlong_initcond/ init_cond

    double precision ring_inner, ring_outer, center(2)
    common /latlong_initcond_parms2/ ring_inner, ring_outer, center

    integer ring_units
    common /latlong_initcond_parms3/ ring_units

    integer mcapa
    common /com_mcapa/ mcapa

    double precision maxelev
    common /com_extruded/ maxelev

    double precision longitude(2), latitude(2), scale(3)

    pi = 4.d0*atan(1.d0)
    pi2 = 2.0*pi

    open(10,file='setprob.data')
    read(10,*) example
    read(10,*) mcapa
    read(10,*) init_cond
    read(10,*) revs_per_second

    read(10,*) longitude(1)
    read(10,*) longitude(2)
    read(10,*) latitude(1)
    read(10,*) latitude(2)

    read(10,*) maxelev

    read(10,*) scale(1)
    read(10,*) scale(2)
    read(10,*) scale(3)


    read(10,*) center(1)
    read(10,*) center(2)

    read(10,*) ring_inner
    read(10,*) ring_outer
    read(10,*) ring_units
    close(10)

    deg2rad = pi/180.d0


end
