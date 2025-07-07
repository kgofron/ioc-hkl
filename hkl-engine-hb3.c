#include <gsl/gsl_sys.h>                // for gsl_isnan
#include "hkl-factory-private.h"        // for autodata_factories_, etc
#include "hkl-pseudoaxis-common-hkl-private.h"  // for hkl_mode_operations, etc
#include "hkl-pseudoaxis-common-psi-private.h"  // for hkl_engine_psi_new, etc
#include "hkl-pseudoaxis-common-q-private.h"  // for hkl_engine_q2_new, etc
#include "hkl-pseudoaxis-common-tth-private.h"  // for hkl_engine_tth2_new, etc
#include "hkl-pseudoaxis-common-readonly-private.h"

/**************/
/* Axes names */
/**************/

#define OMEGA "omega"
#define CHI "chi"
#define PHI "phi"
#define TTH "tth"

/************/
/* Geometry */
/************/

#define HKL_GEOMETRY_HB3_DESCRIPTION				\
	"+ neutron source fix along the :math:`\\vec{x}` direction (1, 0, 0)\n" \
	"+ 3 axes for the sample\n"					\
	"\n"								\
	"  + **" OMEGA "** : rotating around the :math:`\\vec{z}` direction (0, 0, 1)\n" \
	"  + **" CHI "** : rotating around the :math:`-\\vec{y}` direction (1, 0, 0)\n" \
	"  + **" PHI "** : rotating around the :math:`\\vec{x}` direction (1, 0, 0)\n" \
	"\n"								\
	"+ 1 axes for the detector\n"					\
	"\n"								\
	"  + **" TTH "** : rotation around the :math:`-\\vec{y}` direction (0, 0, 1)\n"

static const char* hkl_geometry_hb3_axes[] = {OMEGA, CHI, PHI, TTH};

static HklGeometry *hkl_geometry_new_hb3(const HklFactory *factory)
{
	HklGeometry *self = hkl_geometry_new(factory, &hkl_geometry_operations_defaults);
	HklHolder *h;

	h = hkl_geometry_add_holder(self);
	hkl_holder_add_rotation(h, OMEGA, 0, 0, 1, &hkl_unit_angle_deg);
	hkl_holder_add_rotation(h, CHI, 1, 0, 0, &hkl_unit_angle_deg);
	hkl_holder_add_rotation(h, PHI, 1, 0, 0, &hkl_unit_angle_deg); // don't need

	h = hkl_geometry_add_holder(self);
	hkl_holder_add_rotation(h, TTH, 0, 0, 1, &hkl_unit_angle_deg);

	return self;
}

/*********/
/* Modes */
/*********/

/*****************/
/* mode readonly */
/*****************/

REGISTER_READONLY_INCIDENCE(hkl_engine_template_incidence_new,
			    P99_PROTECT({OMEGA, CHI, PHI}),
			    surface_parameters_y);

REGISTER_READONLY_EMERGENCE(hkl_engine_template_emergence_new,
			    P99_PROTECT({OMEGA, CHI, PHI, TTH}),
			    surface_parameters_y);

/***********/
/* Engines */
/***********/

static HklEngine *hkl_engine_hb3_hkl_new(HklEngineList *engines)
{
	HklEngine *self;
	HklMode *default_mode;

	self = hkl_engine_hkl_new(engines);

	default_mode = bissector();
	hkl_engine_add_mode(self, default_mode);
	hkl_engine_mode_set(self, default_mode);

    // add modes here

	return self;
}

static HklEngine *hkl_engine_hb3_psi_new(HklEngineList *engines)
{
	HklEngine *self;
	HklMode *default_mode;

	self = hkl_engine_psi_new(engines);

	default_mode = psi();
	hkl_engine_add_mode(self, default_mode);
	hkl_engine_mode_set(self, default_mode);

	return self;
}

/***************/
/* Engine list */
/***************/

static HklEngineList *hkl_engine_list_new_HB3(const HklFactory *factory)
{
	HklEngineList *self = hkl_engine_list_new();

	hkl_engine_hb3_hkl_new(self);
	hkl_engine_hb3_psi_new(self);
	hkl_engine_template_incidence_new(self);
	hkl_engine_template_emergence_new(self);

	return self;
}

/* Register the diffractometer into the factory */
REGISTER_DIFFRACTOMETER(hb3, "HB3", HKL_GEOMETRY_HB3_DESCRIPTION);
