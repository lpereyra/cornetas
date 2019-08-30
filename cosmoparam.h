#ifndef COSMOPARAM_H
#define COSMOPARAM_H

#define RHOCRIT 2.77525E11   // Densidad crítica del Universo [Msol h² / Mpc³]
#define GCONS 6.67300E-20    // Constante de Gravitación [km³ / kg / seg²]    
#define Msol 1.9891E30       // Masa del Sol [Kg]                             
#define Kpc 3.08568025E16    // Kiloparsec -> Kilometro                       

struct cosmoparam
{
    double  omegam			;  /* Omega Materia                         */
    double  omegal			;  /* Omega Lambda                          */
    double  omegak			;  /* Omega Curvatura                       */
		double  hparam      ;  /* Parámetro de Hubble adimensional      */
		double  lbox        ;  /* Lado del box [Kpc / h]                */
		double  Mpart       ;  /* Masa de la partícula [10^10 Msol / h] */
		double  redshift		;  /* Redshift                              */
		double  aexp     		;  /*                                       */
		double  Hubble_a    ;  /*                                       */
		double  soft        ;  /* Softening [kpc / h]                   */
		type_int  npart     ;  /* Número de partículas                  */
  	type_int ngrup_sample  ;  /* Número de halos centros    */
#ifndef HALOS_PARTICULAS
    type_int  ngrup     ;  /* Numero de grupos                      */
#endif
} cp;

#endif
