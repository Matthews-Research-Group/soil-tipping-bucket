#include "module_library.h"
#include "../framework/module_creator.h"  // for create_mc

// Include all the header files that define the modules.
#include "canopy_potential.h"
#include "example_module.h"

creator_map BioCroWater::module_library::library_entries =
{
    {"canopy_potential", &create_mc<canopy_potential>},
    {"example_module", &create_mc<example_module>}
};
